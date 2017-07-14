package it.unich.jandom.domains.numerical.octagon

import it.unich.jandom.domains.numerical._

import scala.language.higherKinds
import scala.language.postfixOps
import spire.math.Rational
import breeze.math.Ring

/**
  * Created by fsestini on 7/11/17.
  *
  * Type of abstract octagons (elements of the octagon abstract domain),
  * parameterized by the underlying DBM.
  *
  * Extremely temporary...
  *
  * NOTE: We enforce that the dbm that is used to create the abstract octagon
  *       is *closed*. Hence we do not need to close it every time we perform
  *       some operation, but we need to ensure that it gets closed after any
  *       operation, if we want to construct a new abstract octagon.
  */
case class AbstractOctagon[M[+_, _]](dbm: M[Closed, Double], e: DifferenceBoundMatrix[M]) {
  def dimension: Int = e.nOfVars(dbm)

  def join(other: AbstractOctagon[M])
    (implicit ifield: InfField[Double]): AbstractOctagon[M] =
    new AbstractOctagon(e.dbmUnion(dbm, other.dbm), e)

  def meet(other: AbstractOctagon[M])
    (implicit ifield: InfField[Double]): AbstractOctagon[M] =
    new AbstractOctagon(e.strongClosure(e.dbmIntersection(dbm, other.dbm)), e)

  def forget(vi: VarIndex): AbstractOctagon[M] =
    new AbstractOctagon(e.forget(vi)(dbm), e)

  def top = new AbstractOctagon(e.topDBM[Double], e: DifferenceBoundMatrix[M])
  def bottom = new AbstractOctagon(e.bottomDBM[Double], e: DifferenceBoundMatrix[M])

  def projectInterval(v: VarIndex, closed: M[Closed, Double]): (Double, Double) = {
    if (e.isBottomDBM(closed)) {
      (Double.PositiveInfinity, Double.NegativeInfinity)
    } else {
      val o = v match { case VarIndex(i) => for {
        p1 <- e.get(2 * i - 1, 2 * i)(closed)
        p2 <- e.get(2 * i, 2 * i - 1)(closed)
      } yield (- p1 / 2, p2 / 2) }
      o match {
        case Some(pair) => pair
        case None => throw new IllegalArgumentException()
      }
    }
  }

  private def createInterval(
    low: Array[Double], high: Array[Double], isEmpty: Boolean)
      : BoxDoubleDomain#Property = {
    val dom: BoxDoubleDomain = BoxDoubleDomain(false)
    new dom.Property(low, high, isEmpty)
  }

  def toInterval: BoxDoubleDomain#Property = {
    val closed = e.strongClosure(dbm)
    val l: List[(Double, Double)] =
      (0 until e.nOfVars(dbm))
        .map(i => projectInterval(VarIndex(i), closed)).toList
    val (low, high) = l.unzip
    createInterval(low.toArray, high.toArray, false)
  }

  private def fromFun[A](d: Int, f: (Int, Int) => A): M[DBMState, A] = ???

  private def forSomeVar(
    vars: Seq[VarIndex])(p: VarIndex => Boolean): Option[VarIndex] =
    squash(vars.map(x => if (p(x)) Some(x) else None).toList)

  def fromInterval(box: BoxDoubleDomain#Property): AbstractOctagon[M] = {
    val indices = (0 until dimension).map(x => VarIndex(x))
    val chooser = forSomeVar(indices) _
    val f: (Int, Int) => Double = (i, j) => {
      val g1: VarIndex => Boolean = k => i == 2 * k.i && j == 2 * k.i - 1
      val g2: VarIndex => Boolean = k => j == 2 * k.i && i == 2 * k.i - 1
      (chooser(g1), chooser(g2)) match {
        case (Some(VarIndex(k)), _) => 2 * box.asPair(k)._2
        case (None, Some(VarIndex(k))) => - 2 * box.asPair(k)._1
        case (None, None) => Double.PositiveInfinity
      }
    }
    // TODO not sure if we have to strongly close this...
    val newM: M[DBMState, Double] = fromFun(box.dimension, f)
    new AbstractOctagon(e.strongClosure(newM), e)
  }

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

  def decideLinearForm(assignedVar: VarIndex, lf: LinearForm)
      : Option[ExactLinearForm] =
    lf.pairs.toList match {
      case Nil => Some(ConstExact(lf.known))
      case ((other, coeff) :: Nil) =>
        (assignedVar == other, coeff == -1, coeff == 1) match {
          case (true, true, _) => Some(SingleExact(Negative, lf.known))
          case (true, _, true) => Some(SingleExact(Positive, lf.known))
          case (false, true, _) => Some(DoubleExact(VarIndex(other), Negative, lf.known))
          case (false, _, true) => Some(DoubleExact(VarIndex(other), Positive, lf.known))
          case _ => None
        }
      case _ => None
    }


  def singleConstantExactAssignment(v: VarIndex, const: Double)
    (dbm: M[DBMState, Double]): M[DBMState, Double] = {
    val f: (Int, Int) => Double = (i, j) =>
      if (i == 2 * v.i - 1 && j == 2 * v.i) -2 * const else
        if (i == 2 * v.i && j == 2 * v.i - 1) 2 * const else
          forceOption(e.get(i, j)(e.forget(v)(dbm)))
    e.update(f)(dbm)
  }

  // single exact assignments preserve strong closure
  def singlePositiveExactAssignment[S <: DBMState](v: VarIndex, const: Double)
    (dbm: M[S, Double]) : M[S, Double] = {
    e.addScalarOnVar(v, const)(dbm)
    // CORRECT IMPLEMENTATION:
    // val f: (Int, Int) => Double = (i, j) => {
    //   val g1 = (i == 2 * v - 1 && j != 2 * v - 1 && j != 2 * v) ||
    //            (j == 2 * v && i != 2 * v - 1 && i != 2 * v)
    //   val g2 = (i != 2 * v - 1 && i != 2 * v && j == 2 * v - 1) ||
    //            (j != 2 * v - 1 && j != 2 * v && i == 2 * v)
    //   val g3 = i == 2 * v - 1 && j == 2 * v
    //   val g4 = i == 2 * v && j == 2 * v - 1
    //   if (g1) forceOption(e.get(i, j)(dbm)) - const else
    //     if (g2) forceOption(e.get(i, j)(dbm)) + const else
    //       if (g3) forceOption(e.get(i, j)(dbm)) - 2 * const else
    //         if (g4) forceOption(e.get(i, j)(dbm)) + 2 * const else
    //           forceOption(e.get(i, j)(dbm))
    // }
    // e.update(f)(dbm)
  }

  def doublePositiveExactAssignment(
    vi: VarIndex, vother: VarIndex, const: Double)
    (dbm: M[DBMState, Double]): M[DBMState, Double] = {
    val v = vi.i
    val other = vother.i
    val f: (Int, Int) => Double = (i, j) => {
      val g1 = (i == 2 * v - 1 && j == 2 * other - 1) ||
        (i == 2 * other && j == 2 * v)
      val g2 = (i == 2 * other - 1 && j == 2 * v - 1) ||
        (i == 2 * v && j == 2 * other)
      if (g1) (- const) else
        if (g2) const else
          forceOption(e.get(i, j)(e.forget(vi)(e.strongClosure(dbm))))
    }
    e.update(f)(dbm)
  }

  private def signed(i: Int): Int = ???

  // x := - x
  // this preserves strong closure
  def singleNegativeZeroAssignment[S <: DBMState](v: VarIndex)
    (dbm: M[S, Double]): M[S, Double] = {
    e.flipVar(v)(dbm)
    // CORRECT IMPLEMENTATION:
    // val f: (Int, Int) => Double = (i, j) => {
    //   if (i == 2 * v - 1 || i == 2 * v) {
    //     if (j == 2 * v - 1 || j == 2 * v)
    //       forceOption(e.get(signed(i), signed(i))(dbm))
    //     else
    //       forceOption(e.get(signed(i), j)(dbm))
    //   } else {
    //     if (j == 2 * v - 1 || j == 2 * v)
    //       forceOption(e.get(i, signed(j))(dbm))
    //     else
    //       forceOption(e.get(i, j)(dbm))
    //   }
    // }
    // e.update(f)(dbm)
  }

  // x := - y
  def doubleNegativeZeroAssignment(v: VarIndex, other: VarIndex)
      : M[DBMState, Double] => M[DBMState, Double] =
    doublePositiveExactAssignment(v, other, 0) _ andThen
      singleNegativeZeroAssignment(v)

  // x := - x + c
  def singleNegativeExactAssignment(v: VarIndex, const: Double)
      (dbm: M[Closed, Double]): M[Closed, Double] =
    singlePositiveExactAssignment(v, const)(
      singleNegativeZeroAssignment(v)(dbm))

  // x := - y + c
  def doubleNegativeExactAssignment(v: VarIndex, other: VarIndex, const: Double)
      : M[DBMState, Double] => M[DBMState, Double] =
    doubleNegativeZeroAssignment(v, other) andThen
      singlePositiveExactAssignment(v, const)

  private def varLf(v: VarIndex, dimension: Int): LinearForm = {
    val sss: Seq[Rational] =
      (0 to dimension).map(x => if (x == v.i + 1) Rational(1) else Rational(0))
    new DenseLinearForm(sss)
  }

  private def squash[A](l: List[Option[A]]): Option[A] =
    l match {
      case Nil => None
      case (Some(x) :: xs) => Some(x)
      case (None :: xs) => squash(xs)
    }

  def thruIntervals(vi: VarIndex, lf: LinearForm, dimension: Int)
    (dbm: M[Closed, Double]): M[DBMState, Double] = {
    val v = vi.i
    val f: (Int, Int) => Double = (i, j) => {
      if (i == 2 * v && j == 2 * v - 1) {
        val p = lfAsInterval(vi, lf)
        2 * math.max(p _1, p _2)
      } else if (i == 2 * v - 1 && j == 2 * v) {
        val p = lfAsInterval(vi, lf)
        - 2 * math.max(p _1, p _2)
      } else {
        val g1: VarIndex => Boolean = other =>
        (v != other) && ((i == 2 * other.i - 1 && j == 2 * v - 1) ||
          (i == 2 * v && j == 2 * other.i))
        val g2: VarIndex => Boolean = other =>
        (v != other) && ((i == 2 * other.i && j == 2 * v - 1) ||
          (i == 2 * v && j == 2 * other.i - 1))
        val g3: VarIndex => Boolean = other =>
        (v != other) && ((i == 2 * v - 1 && j == 2 * other.i - 1) ||
          (i == 2 * other.i && j == 2 * v))
        val g4: VarIndex => Boolean = other =>
        (v != other) && ((i == 2 * other.i - 1 && j == 2 * v) ||
          (i == 2 * v - 1 && j == 2 * other.i))

        val chooser = forSomeVar((0 until dimension).map(VarIndex)) _
        val r = (chooser(g1), chooser(g2), chooser(g3), chooser(g4)) match {
          case (Some(other), _, _, _) => {
            val p = lfAsInterval(vi, lf - varLf(other, dimension))
            Some(math.max(p _1, p _2))
          }
          case (_, Some(other), _, _) => {
            val p = lfAsInterval(vi, lf + varLf(other, dimension))
            Some(math.max(p _1, p _2))
          }
          case (_, _, Some(other), _) => {
            val p = lfAsInterval(vi, varLf(other, dimension) - lf)
            Some(math.max(p _1, p _2))
          }
          case (_, _, _, Some(other)) => {
            val p = lfAsInterval(vi, - lf - varLf(other, dimension))
            Some(math.max(p _1, p _2))
          }
          case _ => None
        }
        r match {
          case Some(x) => x
          case None => forceOption(e.get(i,j)(e.strongClosure(dbm)))
        }
      }
    }
    e.update(f)(dbm)
  }

  def assignment(v: VarIndex, lf: LinearForm): AbstractOctagon[M] = {
    val f: M[Closed, Double] => M[Closed, Double] = decideLinearForm(v, lf) match {
      case Some(ConstExact(const)) =>
        singleConstantExactAssignment(v, const.toDouble) _ andThen
          e.incrementalClosure(v)
      case Some(SingleExact(Positive, const)) =>
        singlePositiveExactAssignment(v, const.toDouble)
      case Some(SingleExact(Negative, const)) =>
        singleNegativeExactAssignment(v, const.toDouble)
      case Some(DoubleExact(other, Positive, const)) =>
        doublePositiveExactAssignment(v, other, const.toDouble) _ andThen
          e.incrementalClosure(v)
      case Some(DoubleExact(other, Negative, const)) =>
        doubleNegativeExactAssignment(v, other, const.toDouble) andThen
          e.incrementalClosure[Double](v)
      case None => thruIntervals(v, lf, dimension) _ andThen
        e.incrementalClosure(v)
    }
    new AbstractOctagon(f(dbm), e)
  }

  // Evaluation of linear assignment using interval arithmetics.
  def lfAsInterval(v: VarIndex, lf: LinearForm): (Double, Double) = ???

  trait OctaVarCoeff
  object Positive extends OctaVarCoeff { }
  object Negative extends OctaVarCoeff { }
  trait ExactLinearForm
  case class ConstExact(const: Rational) extends ExactLinearForm
  case class SingleExact(varCoeff: OctaVarCoeff, const: Rational) extends ExactLinearForm
  case class DoubleExact(other: VarIndex, varCoeff: OctaVarCoeff, const: Rational) extends ExactLinearForm

}
