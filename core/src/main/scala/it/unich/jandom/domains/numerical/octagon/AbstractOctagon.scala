package it.unich.jandom.domains.numerical.octagon

import it.unich.jandom.domains.numerical._
import scala.language.higherKinds
import spire.math.Rational

/**
  * Created by fsestini on 7/11/17.
  *
  * Type of abstract octagons (elements of the octagon abstract domain),
  * parameterized by the underlying DBM.
  *
  * Extremely temporary...
  */
case class AbstractOctagon[M[_]](dbm: M[Double], e: DifferenceBoundMatrix[M]) {
  def join(other: AbstractOctagon[M]): AbstractOctagon[M] =
    new AbstractOctagon(e.union(e.strongClosure(dbm), e.strongClosure(other.dbm)), e) // check if troo

  def meet(other: AbstractOctagon[M]): AbstractOctagon[M] =
    new AbstractOctagon(e.intersection(dbm, other.dbm), e) // check if troo

  def forget(): AbstractOctagon[M] = ???
  def forget(vi: VarIndex): AbstractOctagon[M]
  def forgetDBM(vi: VarIndex, m: M[Double]): M[Double] = ???
    // strongly close m, then apply forget operator

  // etc...

  def top = new AbstractOctagon(e.topDBM[Double], e: DifferenceBoundMatrix[M])
  def bottom = new AbstractOctagon(e.bottomDBM[Double], e: DifferenceBoundMatrix[M])

  def toInterval: BoxDoubleDomain#Property
  def fromInterval(box: BoxDoubleDomain#Property): AbstractOctagon[M]

  private def forceOption[A](o: Option[A]): A = o match {
    case Some(x) => x
    case None => throw new IllegalArgumentException()
  }

  //////////////////////// ASSIGNMENTS /////////////////////////////////////////

  // There are only two assignment forms which have an exact abstraction in the
  // octagon domain: (x := c) and (x := +/- y + c).
  // Assignments (x := +/- x + c) do not require strongly closed matrix arguments,
  // but preserve the strong closure.
  // Assignments (x := +/- y + c) with x =/= y require a strongly closed argument
  // due to the embedded forget operator. The result is not strongly closed, but
  // can be strongly closed by merely performing an incremental strong closure
  // with respect to the assigned variable x.
  // Thus, a transfer function that keeps matrices in strongly closed form can be
  // computed in quadratic time, in the worst case.

  def decideLinearForm(assignedVar: Int, lf: LinearForm)
      : Option[ExactLinearForm] =
    lf.pairs.toList match {
      case Nil => Some(ConstExact(lf.known))
      case ((other, coeff) :: Nil) =>
        (assignedVar == other, coeff == -1, coeff == 1) match {
          case (true, true, _) => Some(SingleExact(Negative, lf.known))
          case (true, _, true) => Some(SingleExact(Positive, lf.known))
          case (false, true, _) => Some(DoubleExact(other, Negative, lf.known))
          case (false, _, true) => Some(DoubleExact(other, Positive, lf.known))
          case _ => None
        }
      case _ => None
    }

  def singleConstantExactAssignment(v: VarIndex, const: Double)(dbm: M[Double])
      : M[Double] = {
    val f: (Int, Int) => Double = (i, j) =>
      if (i == 2 * v - 1 && j == 2 * v) (- 2 * const) else
        if (i == 2 * v && j == 2 * v - 1) (2 * const) else
          forceOption(e.get(i, j)(forgetDBM(v, dbm)))
    e.update(f)(dbm)
    // new AbstractOctagon(e.update(f)(dbm), e)
  }

  def singlePositiveExactAssignment(v: VarIndex, const: Double)
    (dbm: M[Double]) : M[Double] = {
    val f: (Int, Int) => Double = (i, j) => {
      val g1 = (i == 2 * v - 1 && j != 2 * v - 1 && j != 2 * v) ||
               (j == 2 * v && i != 2 * v - 1 && i != 2 * v)
      val g2 = (i != 2 * v - 1 && i != 2 * v && j == 2 * v - 1) ||
               (j != 2 * v - 1 && j != 2 * v && i == 2 * v)
      val g3 = i == 2 * v - 1 && j == 2 * v
      val g4 = i == 2 * v && j == 2 * v - 1
      if (g1) (forceOption(e.get(i,j)(dbm)) - const) else
        if (g2) (forceOption(e.get(i,j)(dbm)) + const) else
          if (g3) (forceOption(e.get(i,j)) - 2 * const) else
            if (g4) (forceOption(e.get(i, j)(dbm)) + 2 * const) else
              forceOption(e.get(i, j)(dbm))
    }
    e.update(f)(dbm)
    // new AbstractOctagon(e.update(f)(dbm), e)
  }

  def doublePositiveExactAssignment(v: VarIndex, other: VarIndex, const: Double)
    (dbm: M[Double]): M[Double] = {
    val f: (Int, Int) => Double = (i, j) => {
      val g1 = (i == 2 * v - 1 && j == 2 * other - 1) ||
        (i == 2 * other && j == 2 * v)
      val g2 = (i == 2 * other - 1 && j == 2 * v - 1) ||
        (i == 2 * v && j == 2 * other)
      if (g1) (- const) else
        if (g2) const else
          forceOption(e.get(i, j)(forgetDBM(v, e.strongClosure(dbm))))
    }
    e.update(f)(dbm)
    // new AbstractOctagon(e.update(f)(dbm), e)
  }

  private def signed(i: Int): Int = ???

  // x := - x
  def singleNegativeZeroAssignment(v: VarIndex)(dbm: M[Double]): M[Double] = {
    val f: (Int, Int) => Double = (i, j) => {
      if (i == 2 * v - 1 || i == 2 * v) (
        if (j == 2 * v - 1 || j ==  2 * v)
          e.get(signed(i), signed(i))(dbm)
        else
          e.get(signed(i), j)(dbm)
      )
      else (
        if (j == 2 * v - 1 || j ==  2 * v)
          e.get(i, signed(j))(dbm)
        else
          e.get(i, j)(dbm)
      )
    }
    e.update(f)(dbm)
    // new AbstractOctagon(e.update(f)(dbm), e)
  }

  // x := - y
  def doubleNegativeZeroAssignment(v: VarIndex, other: VarIndex): M[Double] => M[Double] =
    doublePositiveExactAssignment(v, other, 0) _ andThen
      singleNegativeZeroAssignment(v) _

  // x := - x + c
  def singleNegativeExactAssignment(v: Int, const: Double): M[Double] => M[Double] =
    singleNegativeZeroAssignment(v) _ andThen
      singlePositiveExactAssignment(v, const) _

  // x := - y + c
  def doubleNegativeExactAssignment(v: Int, other: VarIndex, const: Double): M[Double] => M[Double] =
    doubleNegativeZeroAssignment(v, other) andThen
      singlePositiveExactAssignment(v, const)

  private def sumVar(v: VarIndex, coeff: Rational)(lf: LinearForm): LinearForm = {
    val sss: Seq[Rational] =
      (0 to lf.dimension).map(x => if (x == v + 1) coeff else Rational(0))
    val vlf: LinearForm = new DenseLinearForm(sss)
    lf + vlf
  }

  private def inverseLf(lf: LinearForm): LinearForm = lf * (-1)

  def thruIntervals(v: VarIndex, lf: LinearForm)(dbm: M[Double]): M[Double] = {
    val f: (Int, Int) => Double = (i, j) => {
      val g1 = i == 2 * v && j == 2 * v - 1
      val g2 = i == 2 * v - 1 && j == 2 * v
      // val g3 = i == 2 * v - 1 && j == 2 * v - 1 && i
      ???
    }
    ???
  }

  def assignment(v: VarIndex, lf: LinearForm): AbstractOctagon[M] = {
    val f: M[Double] => M[Double] = decideLinearForm(v, lf) match {
      case Some(ConstExact(const)) =>
        singleConstantExactAssignment(v, const.toDouble)
      case Some(SingleExact(Positive, const)) =>
        singlePositiveExactAssignment(v, const.toDouble)
      case Some(SingleExact(Negative, const)) =>
        singleNegativeExactAssignment(v, const.toDouble)
      case Some(DoubleExact(other, Positive, const)) =>
        doublePositiveExactAssignment(v, other, const.toDouble) _ andThen
          e.incrementalClosure(v)
      case Some(DoubleExact(other, Negative, const)) =>
        doubleNegativeExactAssignment(v, other, const.toDouble) andThen
          e.incrementalClosure(v)
      case None => ???
    }
    new AbstractOctagon(f(dbm), e)
  }

  // Evaluation of linear assignment using interval arithmetics.
  def lfAsInterval(v: VarIndex, lf: LinearForm): (Double, Double) = ???

  type VarIndex = Int

  trait OctaVarCoeff
  object Positive extends OctaVarCoeff { }
  object Negative extends OctaVarCoeff { }
  trait ExactLinearForm
  case class ConstExact(const: Rational) extends ExactLinearForm
  case class SingleExact(varCoeff: OctaVarCoeff, const: Rational) extends ExactLinearForm
  case class DoubleExact(other: VarIndex, varCoeff: OctaVarCoeff, const: Rational) extends ExactLinearForm

}
