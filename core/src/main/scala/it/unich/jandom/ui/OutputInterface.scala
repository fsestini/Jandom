/**
 * Copyright 2013, 2016 Jandom Team
 *
 * This file is part of JANDOM: JVM-based Analyzer for Numerical DOMains
 * JANDOM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * JANDOM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of a
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with JANDOM.  If not, see <http://www.gnu.org/licenses/>.
 */

package it.unich.jandom.ui

import java.io.FileInputStream
import java.io.IOException
import java.nio.file._
import java.nio.file.attribute.BasicFileAttributes

import scala.collection.JavaConversions._
import scala.util.Try

import jdk.internal.org.objectweb.asm.ClassReader
import jdk.internal.org.objectweb.asm.tree.ClassNode
import jdk.internal.org.objectweb.asm.tree.MethodNode

import it.unich.jandom.domains.numerical.NumericalDomain
import it.unich.jandom.domains.objects.ObjectDomainFactory
import it.unich.jandom.targets.parameters.WideningSpecs._
import it.unich.jandom.targets.parameters.NarrowingSpecs._
import it.unich.jandom.parsers.FastParser
import it.unich.jandom.parsers.RandomParser
import it.unich.jandom.targets._
import it.unich.jandom.targets.jvmasm._
import it.unich.jandom.targets.jvmsoot._
import it.unich.jandom.targets.lts._
import it.unich.jandom.targets.slil._
import soot.Scene
import soot.toolkits.graph.Block
import it.unich.scalafix.Box

/**
 * An output interface is a collection of methods for implementing an external interface.
 * @author Gianluca Amato <gianluca.amato@unich.it>
 * @author Francesca Scozzari <francesca.scozzari@unich.it>
 */
object OutputInterface {

  def getWideningStrategies = {
    WideningScopes.values.map(x => x.name)
  }

  def getWideningStrategiesTips = {
    WideningScopes.values.map(x => x.description)
  }

  def getNarrowingStrategies = {
    NarrowingStrategies.values.map(x => x.name)
  }

  def getNarrowingStrategiesTips = {
    NarrowingStrategies.values.map(x => x.description)
  }

  def getNumericalDomains = {
    NumericalDomains.values.map(x => x.name)
  }

  def getNumericalDomainsTips = {
    NumericalDomains.values.map(x => x.description)
  }

  def getObjectDomains = {
    ObjectDomains.values.map(x => x.name)
  }
  def getObjectDomainsTips = {
    ObjectDomains.values.map(x => x.description)
  }

  def getDebugTip = "Generate debug information - for developers only"
  def getWideningDelayTip = "Apply the widening with delay"
  def getRadioBafTip = "Analysis of Java bytecode thorugh the Baf representation of the Soot library"
  def getRadioJimpleTip = "Analysis of Java bytecode thorugh the Jimple representation of the Soot library"
  def getRadioNumericalTip = "Analysis of numerical properties"
  def getRadioObjectTip = "Analysis of object properties"
  def getClassPathTip = "Choose the classpath of the program to analyze"
  def getClassTip = "Choose the class to analyze"
  def getMethodTip = "Choose the method to analyze"
  def getIRTypeTip = "Type of intermediate representation to use"
  def getAnalysisTypeTip = "Choose between an analysis of numerical properties or object-related properties"

  def setParameters[T <: Target[T]](params: Parameters[T], wideningIndex:Int, narrowingIindex:Int, delay:Int, debug:Boolean) {
    params.wideningScope = WideningScopes.values(wideningIndex).value
    params.narrowingStrategy = NarrowingStrategies.values(narrowingIindex).value
    params.widening = DelayedWidening(DefaultWidening, delay)
    params.narrowing = DelayedNarrowing(TrivialNarrowing, 2)
    if (debug) params.debugWriter = new java.io.StringWriter
  }

  /**
   * Analyze a class using Soot.
   * @param method the method to be analyzed
   * @param domain the abstract domain to be used (either numerical or object)
   * @param wideningIndex the widening strategy
   * @param narrowingIndex the narrowing strategy
   * @param delay the widening delay
   * @param debug is true when the debug is active
   * @return a string with the program annotated with the analysis result
   */
  private def analyze[T <: SootCFG[T, Block]](method: SootCFG[T, Block], domain: Any, wideningIndex: Int,
    narrowingIndex: Int, delay: Int, debug: Boolean): String = {
    try {
      val sootScene = Scene.v()
      sootScene.loadBasicClasses()
      val om = new SootObjectModel(sootScene)
      val sootDomain: SootFrameDomain = domain match {
        case domain: NumericalDomain => new SootFrameNumericalDomain(domain)
        case domain: ObjectDomainFactory => new SootFrameObjectDomain(domain(om))
      }
      val tMethod = method.asInstanceOf[T]
      val params = new Parameters[T] { val domain = sootDomain }
      setParameters(params, wideningIndex, narrowingIndex, delay, debug)
      val inte = new TopSootInterpretation[T, params.type](params)
      params.interpretation = Some(inte)
      val ann = tMethod.analyze(params)
      tMethod.mkString(params)(ann)
    } catch {
      case e: UnsupportedSootUnitException =>
        e.getMessage + " : " + e.unit + " Error in analysing bytecode"
      case e: Exception =>
        e.getMessage + " Error in parsing source code"
    }
  }

  def analyze(dir: Path, klass: Int, method: Int, isNumerical: Boolean, isBaf: Boolean, domain: Int, wideningIndex: Int,
    narrowingIndex: Int, delay: Int, debug: Boolean): String = {
    val methods = getSootMethods(dir, klass)
    val selectedMethod = methods.get(method)

    val aDomain = if (isNumerical)
      NumericalDomains.values(domain).value
    else
      ObjectDomains.values(domain).value
    if (isBaf)
      analyze(new BafMethod(selectedMethod), aDomain, wideningIndex, narrowingIndex, delay, debug)
    else
      analyze(new JimpleMethod(selectedMethod), aDomain, wideningIndex, narrowingIndex, delay, debug)
  }

  private def getScene(dir: Path) = {
    val scene = Scene.v()
    scene.loadBasicClasses()
    scene.setSootClassPath(scene.defaultClassPath + java.io.File.pathSeparator + dir.toString)
    scene
  }

  def unzipJar(dir: String, jarFile: String, destDir: String) = {
    val jar = new java.util.jar.JarFile(dir + java.io.File.separator + jarFile)
    val enum = jar.entries()
    while (enum.hasMoreElements()) {
      val file = enum.nextElement(); //(java.util.jar.JarEntry)
      val f = new java.io.File(destDir + java.io.File.separator + file.getName())
      if (file.isDirectory()) { // if its a directory, create it
        f.mkdir()
      } else {
        val is = jar.getInputStream(file); // get the input stream
        val fos = new java.io.FileOutputStream(f);
        while (is.available() > 0) { // write contents of 'is' to 'fos'
          fos.write(is.read())
        }
        fos.close()
        is.close()
      }
    }
  }

  def getMethods(dir: Path, klassIndex: Int) = {
    val scene = getScene(dir)
    val sootKlass = scene.loadClassAndSupport(getClasses(dir)(klassIndex))
    sootKlass.getMethods().map(x => x.getName())
  }

  private def getSootMethods(dir: Path, klassIndex: Int) = {
    val scene = getScene(dir)
    val sootKlass = scene.loadClassAndSupport(getClasses(dir)(klassIndex))
    sootKlass.setApplicationClass()
    sootKlass.getMethods()
  }

  def getClasses(dir: Path): Seq[String] = {
    if (Files.isDirectory(dir)) {
      val fileProcessor = new ClassFileVisitor(dir)
      if (Try(Files.walkFileTree(dir, fileProcessor)).isSuccess) {
        // these two lines are a mess because Scala Swing does not play well with Java 1.7
        fileProcessor.classNameList
      } else Seq[String]()
    } else {
      val classNames = scala.collection.mutable.SortedSet[String]()
      val jar = new java.util.jar.JarFile(dir.toFile)
      val enum = jar.entries()
      while (enum.hasMoreElements()) {
        val file = enum.nextElement().getName
        if (file.endsWith(".class")) classNames += file stripSuffix ".class"
      }
      classNames.toSeq
    }
  }

  def getSootAbstraction(dir: Path, klassIndex: Int, methodIndex: Int, isBaf: Boolean) = {
    val myMethod = getSootMethods(dir, klassIndex).get(methodIndex)
    if (isBaf)
      new BafMethod(myMethod).toString
    else
      new JimpleMethod(myMethod).toString
  }

  def getASMAbstraction(dir: String, klassName: String, methodIndex: Int) = {
    new AsmMethod(getASMMethodsList(dir, klassName).get(methodIndex)).toString
  }

  private def getASMMethodsList(dir: String, klassName: String) = {
    val file = dir + java.io.File.separator + klassName
    val is = new FileInputStream(file)
    val cr = new ClassReader(is)
    val node = new ClassNode()
    cr.accept(node, ClassReader.SKIP_DEBUG)
    val methodList = node.methods.asInstanceOf[java.util.List[MethodNode]]
    methodList
  }

  def getASMMethods(dir: String, klassName: String) = {
    getASMMethodsList(dir, klassName) map { _.name }
  }

  def analyzeASM(dir: String, klassName: String, methodIndex: Int, domain: Int, wideningIndex: Int,
    narrowingIndex: Int, delay: Int, debug: Boolean) = {
    try {
      val numericalDomain = NumericalDomains.values(domain).value
      val jvmDomain = new JVMEnvFixedFrameDomain(numericalDomain)
      val method = new AsmMethod(getASMMethodsList(dir, klassName).get(methodIndex))
      val params = new Parameters[AsmMethod] { val domain = jvmDomain }
      setParameters(params, wideningIndex, narrowingIndex, delay, debug)
      val ann = method.analyze(params)
      method.mkString(ann)
    } catch {
      case e: UnsupportedASMInsnException =>
        e.getMessage + " : " + e.node + " Error in analysing bytecode"
      case e: Exception =>
        e.getMessage + " Error"
    }
  }

  def getRandomText(dir: String, file: String) = {
    val filename = new java.io.File(dir + java.io.File.separator + file);
    try {
      scala.io.Source.fromFile(filename).mkString
    } catch {
      case e: IOException => "File not found"
    }
  }

  def analyzeRandomStr(program: String, domain: Int, widening: Int, narrowing: Int, delay: Int, debug: Boolean) = {
    val parser = RandomParser()
    val numericalDomain = NumericalDomains.values(domain).value
    val result = parser.parseProgram(program) match {
      case parser.Success(program, _) =>
        val params = new Parameters[SLILTarget] { val domain = numericalDomain }
        setParameters(params, widening, narrowing, delay, debug)
        val ann = program.analyze(params)
        params.debugWriter.toString + program.mkString(ann)
      case parser.NoSuccess(msg, next) =>
        msg + " in line " + next.pos.line + " column " + next.pos.column + " Error in parsing source code"
      case _ => "Error in parsing"
    }
    result
  }

  def analyzeRandom(dir: String, file: String, domain: Int, widening: Int, narrowing: Int, delay: Int, debug: Boolean) = {
    try {
      val program = getRandomText(dir, file)
      analyzeRandomStr(program, domain, widening, narrowing, delay, debug)
    } catch {
      case e: IOException => "I/O error"
    }
  }

  def analyzeFastModelStr(program: String, domain: Int, widening: Int, narrowing: Int, delay: Int, debug: Boolean) = {
    try {
      val parser = FastParser()
      val numericalDomain = NumericalDomains.values(domain).value

      val result = parser.parse(program) match {
        case parser.Success(program, _) =>
          val params = new Parameters[LTS] { val domain = numericalDomain }
          setParameters(params, widening, narrowing, delay, debug)
          val ann = program.analyze(params)
          val grafo = program.toDot
          (grafo, params.debugWriter.toString + program.mkString(ann) )
        case parser.NoSuccess(msg, next) =>
          msg + " in line " + next.pos.line + " column " + next.pos.column + " Error in parsing source code"
        case _ => "Error in parsing"
      }
      result
    } catch {
      case e: IOException => "I/O error"
    }
  }
}

private class ClassFileVisitor(rootPath: Path) extends SimpleFileVisitor[Path] {
  private val privateClassNamesList = scala.collection.mutable.SortedSet[String]()
  def classNameList = privateClassNamesList.toSeq
  override def visitFile(aFile: Path, aAttrs: BasicFileAttributes): FileVisitResult = {
    val relativePath = rootPath.relativize(aFile)
    val className = (relativePath.head.toString /: relativePath.tail)(_ + "." + _.toString)
    if (className endsWith ".class")
      privateClassNamesList += className stripSuffix ".class"
    else if (className endsWith ".java")
      privateClassNamesList += className stripSuffix ".java"
    FileVisitResult.CONTINUE;
  }
}
