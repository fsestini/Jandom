//*** Libraries

libraryDependencies ++= Seq(
  "org.apache.commons" % "commons-lang3" % "3.5",
  "org.scalatest" %% "scalatest" % "3.0.2" % "test",
  "org.scalacheck" %% "scalacheck" % "1.13.3" % "test",
  "org.mockito" % "mockito-core" % "2.2.9" % "test",
  "org.spire-math" %% "spire" % "0.12.0",
  "it.unich.scalafix" %% "scalafix" % "0.6.0",
  "org.rogach" %% "scallop" % "2.0.3",
  "org.scala-lang.modules" %% "scala-swing" % "1.0.2",
  "org.scala-lang.modules" %% "scala-parser-combinators" % "1.0.4",
  "org.scala-lang" % "scala-reflect" % scalaVersion.value,
  // ASM is included in the Soot Jar
  "ca.mcgill.sable" % "soot" %"3.0.0-SNAPSHOT",
  "org.scalaz" %% "scalaz-core" % "7.2.14",
  "org.scala-graph" %% "graph-core" % "1.11.5",
  "org.scalanlp" %% "breeze" % "0.12",
  "org.scalanlp" %% "breeze-viz" % "0.12"
)

updateOptions := updateOptions.value.withLatestSnapshots(false)

//*** Additional source directories for PPL

unmanagedJars in Compile ++= (pplJar.value map file).toSeq

unmanagedSourceDirectories in Compile ++= (pplJar.value map { _ => (sourceDirectory in Compile).value / "ppl" }).toSeq

unmanagedSourceDirectories in Test ++= (pplJar.value map { _ => (sourceDirectory in Test).value / "ppl" }).toSeq

//*** IntelliJ Idea

ideOutputDirectory in Compile := Some(new File("core/target/idea/classes"))

ideOutputDirectory in Test := Some(new File("core/target/idea/test-classes"))

//*** Eclipse plugin

EclipseKeys.createSrc := EclipseCreateSrc.Default +  EclipseCreateSrc.Managed

EclipseKeys.executionEnvironment := Some(EclipseExecutionEnvironment.JavaSE17)

EclipseKeys.eclipseOutput := Some("target.eclipse")

// It would be nice to be able to exclude resource directories from compilation.

managedSourceDirectories in Test := Seq()

managedResourceDirectories in Test := Seq()

managedResourceDirectories in Compile := Seq()

unmanagedResourceDirectories in Compile := Seq()

//*** BuildInfo plugin

gitHeadCommitSHA := Process("git rev-parse HEAD").lines.head

buildInfoKeys := Seq[BuildInfoKey](name, version, scalaVersion, sbtVersion, gitHeadCommitSHA)

buildInfoPackage := "it.unich.jandom"


