package it.unich.jandom.domains.numerical.octagon

import breeze.math.Field
import it.unich.jandom.domains.numerical.LinearForm

// Distinguish integers used as variable indices
case class VarIndex(i: Int) extends Ordered[VarIndex] {
  def compare(that: VarIndex) = this.i compare that.i
}

object VarIndexOps {
  sealed trait OctaVarCoeff
  object Positive extends OctaVarCoeff { }
  object Negative extends OctaVarCoeff { }

  def vcAsNumeral[A](c: OctaVarCoeff)(implicit ifield: Field[A]): A = c match {
    case Positive => ifield.one
    case Negative => ifield.inverse(ifield.one)
  }

  def varPlus(v: VarIndex): Int = 2 * v.i
  def varMinus(v: VarIndex): Int = 2 * v.i + 1
  def signed(i: Int): Int = if (i % 2 == 0) i + 1 else i - 1

  def inverseVarPlusMinus(coeff: OctaVarCoeff, dimension: Int, i: Int): Option[VarIndex] =
    if (0 <= i && i < dimension && toIndexAndCoeff(i)._2 == coeff)
      Some(toIndexAndCoeff(i)._1) else None
  def toIndexAndCoeff(i: Int): (VarIndex, OctaVarCoeff) =
    if (i % 2 == 0) (VarIndex(i / 2), Positive)
    else (VarIndex((i - 1) / 2), Negative)
  def fromIndexAndCoeff(vi: VarIndex, c: OctaVarCoeff): Int = c match {
    case Positive => varPlus(vi)
    case Negative => varMinus(vi)
  }
}

object VarIndexUtils {
  def forSomeVar(
                  vars: Seq[VarIndex])(p: VarIndex => Boolean): Option[VarIndex] =
    (vars.map(x => if (p(x)) Some(x) else None).toList).flatten.headOption
  // Evaluation of linear assignment using interval arithmetics.
  def lfAsInterval(v: VarIndex, lf: LinearForm): (Double, Double) = ???
}

// Distinguish integers used as variable count or matrix dimension
case class VarCount(count: Int) { require(count >= 0) }
case class Dimension(dim: Int) { require(dim % 2 == 0 && dim >= 0) }

object CountOps {
  def halvedDimension(d: Dimension): VarCount = VarCount(d.dim / 2)
  def doubledVarCount(c: VarCount): Dimension = Dimension(c.count * 2)
  def addOne(c: VarCount): VarCount = VarCount(c.count + 1)
  def subOne(c: VarCount): VarCount = VarCount(c.count - 1)
  def inDimension(i: Int, j: Int, dim: Dimension): Boolean =
    0 <= i && i < dim.dim && 0 <= j && j < dim.dim
  def grid(d: Dimension): Seq[(Int, Int)] = for {
    i <- 0 until d.dim
    j <- 0 until d.dim
  } yield (i, j)
  def allIndices(d: Dimension): Seq[Int] = 0 until d.dim
  def allVars(vc: VarCount): Seq[VarIndex] = (0 until vc.count).map(VarIndex)
}
