package it.unich.jandom.domains.numerical.octagon

import it.unich.jandom.domains.numerical._

import scala.language.higherKinds
import scala.language.postfixOps
import VarIndexOps._
import spire.math.Rational
import breeze.math.Ring
import it.unich.jandom.utils.numberext.RationalExt
import it.unich.jandom.domains.numerical.octagon.OctagonalConstraint._

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

case class AbstractOctagon[M[_, _]](dbm: M[Closed, Double], e: DifferenceBoundMatrix[M])
  extends NumericalProperty[AbstractOctagon[M]] {

  def dimension: Int = e.nOfVars(dbm)

  def union(other: AbstractOctagon[M]): AbstractOctagon[M] =
    AbstractOctagon(e.dbmUnion(dbm, other.dbm)(InfField.infFieldDouble), e)

  def intersection(other: AbstractOctagon[M]): AbstractOctagon[M] =
    AbstractOctagon(e.strongClosure(e.dbmIntersection(dbm, other.dbm).elem), e)

  def forget(vi: VarIndex): AbstractOctagon[M] =
    AbstractOctagon(e.forget(vi)(dbm), e)

  def top = AbstractOctagon(e.topDBM[Double](e.nOfVars(dbm)), e: DifferenceBoundMatrix[M])
  def bottom = AbstractOctagon(e.bottomDBM[Double](e.nOfVars(dbm)), e: DifferenceBoundMatrix[M])

  def widening(other: AbstractOctagon[M]): AbstractOctagon[M] =
    AbstractOctagon(e.strongClosure(e.widening(dbm, other.dbm).elem), e)

  def narrowing(other: AbstractOctagon[M]): AbstractOctagon[M] =
    AbstractOctagon(e.strongClosure(e.narrowing(dbm, other.dbm).elem), e)

  def projectInterval(v: VarIndex, closed: M[Closed, Double]): (Double, Double) = {
    val maybeInterval: Option[(Double, Double)] = for {
      p1 <- e.get(varPlus(v), varMinus(v))(closed)
      p2 <- e.get(varMinus(v), varPlus(v))(closed)
    } yield (- p1 / 2, p2 / 2)
    maybeInterval match {
      case Some(interval) => interval
      case None => (Double.PositiveInfinity, Double.NegativeInfinity)
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
    createInterval(low.toArray, high.toArray, isEmpty = false)
  }

  private def fromFun[A](d: Int, f: (Int, Int) => A)
  : ExistsDBM[({ type T[S] = M[S, Double]})#T] = ???

  private def forSomeVar(
    vars: Seq[VarIndex])(p: VarIndex => Boolean): Option[VarIndex] =
    squash(vars.map(x => if (p(x)) Some(x) else None).toList)

  def fromInterval(box: BoxDoubleDomain#Property): AbstractOctagon[M] = {
    val indices = (0 until dimension).map(x => VarIndex(x))
    val chooser = forSomeVar(indices) _
    val f: (Int, Int) => Double = (i, j) => {
      val g1: VarIndex => Boolean = k => i == varMinus(k) && j == varPlus(k)
      val g2: VarIndex => Boolean = k => j == varMinus(k) && i == varPlus(k)
      (chooser(g1), chooser(g2)) match {
        case (Some(VarIndex(k)), _) => 2 * box.asPair(k)._2
        case (None, Some(VarIndex(k))) => - 2 * box.asPair(k)._1
        case (None, None) => Double.PositiveInfinity
      }
    }
    // TODO not sure if we have to strongly close this...
    AbstractOctagon(e.strongClosure(fromFun(box.dimension, f).elem), e)
  }

  private def forceOption[A](o: Option[A]): A = o match {
    case Some(x) => x
    case None => throw new IllegalArgumentException()
  }

  //////////////////////// BOOL TESTS //////////////////////////////////////////

  /**
   * Computes intersection with `lf != 0`.
   *
   * Amounts to identity as Mine06 does not give a better abstraction
   * than id
   */
  def linearDisequality (lf: LinearForm) = this

  /**
   * Intersection with the half-plane `lf <= 0`.
   *
   * Computed according to Mine06; particularly, 5 cases for an exact
   * abstraction are given, otherwise the intersection is computed by
   * fallback on intervals.
   *
   * The cases given in Mine06 fig 20 are as follows:
   *
   *  1. ({ V_l + c <= 0 ?}(m))_ij = min(m_ij, -2c) if i=2l, j=2l-1
   *  2. ({-V_l + c <= 0 ?}(m))_ij = min(m_ij, -2c) if i=2l-1, j=2l
   *  3. ({V_l - V_k + c  <= 0 ?}(m))_ij = min(m_ij, -c) if i=2k-1, j=2l-1
   *                                                     or i=2l, j=2k
   *  4. ({V_l + V_k + c  <= 0 ?}(m))_ij = min(m_ij, -c) if i = 2k, j = 2l - 1
   *                                                     or i = 2l, j = 2k - 1
   *  5. ({-V_l - V_k + c <= 0 ?}(m))_ij = min(m_ij, -c) if i = 2k - 1, j = 2l
   *                                                     or i = 2l - 1, j = 2k
   *  6. Fallback on interval-based abstraction
   *
   * Notes:
   *  - [a,b] in Mine06 replaced by c, we don't do range assignments
   *  - j0, i0 replaced by k, l for readability
   *  - Always m_ij if not specified
   */
  def linearInequalityEx (lf: LinearForm): ExistsMDouble = {
    sealed abstract class AbstractTest
    case class Fallback() extends AbstractTest
    sealed abstract class ExactTest extends AbstractTest
    case class Case1Test(val vl : Int, val c : Rational) extends ExactTest
    case class Case2Test(val vl : Int, val c : Rational) extends ExactTest
    case class Case3Test(val vl : Int, val vk : Int, val c : Rational) extends ExactTest
    case class Case4Test(val vl : Int, val vk : Int, val c : Rational) extends ExactTest
    case class Case5Test(val vl : Int, val vk : Int, val c : Rational) extends ExactTest
    /**
      * Decodes a Linearform into an AbstractTest, ie one of 6 cases discussed in Mine06.
      */
    def decodeTest (lf : LinearForm) : AbstractTest =
      if ( !lf.homcoeffs.exists { !List(1,-1).contains(_) }
        | lf.homcoeffs.size > 2) {
        Fallback()
        // Mine06 fig 20 does not give an exact abstraction for > 2
        // coeffs or coeffs != 1, -1
      } else if (lf.coeffs.size == 0) {
        throw new IllegalArgumentException("???")
        // Surely 0+0+...+0<=0?
        // TODO: perhaps this is a legit use case, consider removing
      } else {
        val c = lf.known
        if (lf.homcoeffs.size == 1) {
          // Cases 1, 2
          val vl = lf.pairs.head
          if (vl._2 == 1)
            Case1Test (vl._1, c)
          else {
            assert(vl._2 == -1)
            Case2Test (vl._1, c)
          }
        } else if (lf.homcoeffs.size == 2) {
          if (lf.homcoeffs.exists(_ == 1)) {
            // Cases 3, 4
            val vl = lf.pairs.filter(_._2 == 1).head
            val vk = lf.pairs.diff(List(vl)).head
            if (vk._2 == -1)
              Case3Test (vl._1, vk._1, c)
            else
              Case4Test (vl._1, vk._1, c)
          } else {
            // Case 5
            val vl = lf.pairs.head
            val vk = lf.pairs.diff(List(vl)).head
            Case5Test (vl._1, vk._1, c)
          }
        } else {
          // TODO: temporary, eventually remove
          throw new RuntimeException("If you got here my logic is broken, sorry...")
        }
      }

    val t : AbstractTest = decodeTest(lf)
    t match {
      case Fallback() => {
        val lh = fromInterval(toInterval.linearInequality(lf))
        def g (i : Int, j : Int) : Double =
          if (j == i-1 & i%2 == 0 & i/2 <= lf.dimension) { // Case 1 with i/2 = j0
            val interval = this.toInterval
            2*(interval.maximize(LinearForm.v(i/2) - lf).doubleValue())
          } else if (i == j-1 & j%2 == 0 & j/2 <= lf.dimension) { // Case 2 with j/2 = j0
            val interval = this.toInterval
            2*(interval.maximize(LinearForm.v(j/2) - lf).doubleValue())
          } else if ((i%2 == 1 & j%2 == 1 & i != j & (i+1)/2 <= lf.dimension & (j+1)/2 <= lf.dimension)
                     // Case 3a, (i+1)/2 = i0, (j+1)/2 = j0
                    |(i%2 == 0 & j%2 == 0 & i != j & i/2 <= lf.dimension & j/2 < lf.dimension)
                     // Case 3b, i/2 = i0, j/2 = j0
                    ) {
            val interval = this.toInterval
            val j0 = (j+1)/2
            val i0 = (i+1)/2
            2*(interval.maximize(LinearForm.v(j0) - LinearForm.v(i0) - lf).doubleValue())
          } else if (i%2 == 0 & j%2 == 1 & j != i-1 & (i/2) <= lf.dimension & (j+1)/2 <= lf.dimension) {
            // Case 4, i0 = i/2, j0 = j+1/2
            val interval = this.toInterval
            val j0 = (j+1)/2
            val i0 = i/2
            2*(interval.maximize(LinearForm.v(j0) + LinearForm.v(i0) - lf).doubleValue())
          } else if (i%2 == 1 & j%2 == 0 & i != j-1  & (j/2) <= lf.dimension & (i+1)/2 <= lf.dimension) {
            // Case 5, j0 = j, i0 = (i+1)/2
            val interval = this.toInterval
            val j0 = j/2
            val i0 = (i+1)/2
            2*(interval.maximize(- LinearForm.v(j0) - LinearForm.v(i0) - lf).doubleValue())
          } else {
            // id case
            e.get(i,j)(dbm).get
          }
        val f = (i : Int, j : Int) => math.min(
          e.get(i,j)(dbm).get,
          g(i,j)
        )
        e.update(f)(dbm)
      }
      case et : ExactTest => {
        val f = et match {
          case Case1Test(vl, c) => {
            (i : Int, j : Int) =>
            if (i == 2*vl & j == 2*vl - 1)
              math.min(-2*c.toDouble, e.get(i,j)(dbm).get)
            else
              e.get(i,j)(dbm).get
          }
          case Case2Test(vl, c) => {
            (i : Int, j : Int) =>
            if (i == 2*vl-1& j == 2*vl)
              math.min(-2*c.toDouble, e.get(i,j)(dbm).get)
            else
              e.get(i,j)(dbm).get
          }
          case Case3Test(vl, vk, c) => {
            (i : Int, j : Int) =>
            if ((i == 2*vk-1 & j == 2*vl-1)
              | (i == 2*vk & j == 2*vl)
            )
              math.min(-1*c.toDouble, e.get(i,j)(dbm).get)
            else
              e.get(i,j)(dbm).get
          }
          case Case4Test(vl, vk, c) => {
            (i : Int, j : Int) =>
            if ((i == 2*vk & j == 2*vl-1)
              | (i == 2*vl & j == 2*vk-1)
            )
              math.min(-1*c.toDouble, e.get(i,j)(dbm).get)
            else
              e.get(i,j)(dbm).get
          }
          case Case5Test(vl, vk, c) => {
            (i : Int, j : Int) =>
            if ((i == 2*vk-1 & j == 2*vl)
              | (i == 2*vl-1 & j == 2*vk)
            )
              math.min(-1*c.toDouble, e.get(i,j)(dbm).get)
            else
              e.get(i,j)(dbm).get
          }
        }
        e.update(f)(dbm)
      }
    }
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
        (assignedVar.i == other, coeff == -1, coeff == 1) match {
          case (true, true, _) => Some(SingleExact(Negative, lf.known))
          case (true, _, true) => Some(SingleExact(Positive, lf.known))
          case (false, true, _) => Some(DoubleExact(VarIndex(other), Negative, lf.known))
          case (false, _, true) => Some(DoubleExact(VarIndex(other), Positive, lf.known))
          case _ => None
        }
      case _ => None
    }

  type ExistsMDouble = ExistsDBM[({ type T[S] = M[S, Double]})#T]

  def singleConstantExactAssignment[S <: DBMState]
    (v: VarIndex, const: Double)
    (dbm: M[S, Double]): ExistsMDouble = {

    val f: (Int, Int) => Double = (i, j) =>
      if (i == varPlus(v) && j == varMinus(v)) -2 * const else
        if (i == varMinus(v) && j == varPlus(v)) 2 * const else
          forceOption(e.get(i, j)(e.forget(v)(dbm)))
    e.update(f)(dbm)
  }

  // single exact assignments preserve strong closure
  def singlePositiveExactAssignment[S <: DBMState]
    (v: VarIndex, const: Double)
    (dbm: M[S, Double]) : M[S, Double] = e.addScalarOnVar(v, const)(dbm)

  def doublePositiveExactAssignment[S <: DBMState]
    (vi: VarIndex, vother: VarIndex, const: Double)
    (dbm: M[S, Double]): ExistsMDouble = {

    val f: (Int, Int) => Double = (i, j) => {
      val g1 = (i == varPlus(vi) && j == varPlus(vother)) ||
               (i == varMinus(vother) && j == varMinus(vi))
      val g2 = (i == varPlus(vother) && j == varPlus(vi)) ||
               (i == varMinus(vi) && j == varMinus(vother))
      if (g1) -const else
        if (g2) const else
          forceOption(e.get(i, j)(e.forget(vi)(e.strongClosure(dbm))))
    }
    e.update(f)(dbm)
  }

  // x := - x
  // this preserves strong closure
  def singleNegativeZeroAssignment[S <: DBMState]
    (v: VarIndex)(dbm: M[S, Double]): M[S, Double] = e.flipVar(v)(dbm)

  type ExistsM[S] = M[S, Double]

  // x := - y
  def doubleNegativeZeroAssignment[S <: DBMState]
    (v: VarIndex, other: VarIndex)(dbm: M[S, Double]): ExistsMDouble = {
    val m = doublePositiveExactAssignment(v, other, 0)(dbm)
    val mm = singleNegativeZeroAssignment(v)(m.elem)
    MkEx[m.State, ExistsM](mm)
  }

  // x := - x + c
  def singleNegativeExactAssignment(v: VarIndex, const: Double)
                                   (dbm: M[Closed, Double]): M[Closed, Double] =
    singlePositiveExactAssignment(v, const)(
      singleNegativeZeroAssignment(v)(dbm))

  // x := - y + c
  def doubleNegativeExactAssignment[S <: DBMState]
    (v: VarIndex, other: VarIndex, const: Double)
    (dbm: M[S, Double]): ExistsMDouble = {
    val m = doubleNegativeZeroAssignment(v, other)(dbm)
    val mm = singlePositiveExactAssignment(v, const)(m.elem)
    MkEx[m.State, ExistsM](mm)
  }

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
    (dbm: M[Closed, Double]): ExistsMDouble = {
    // val v = vi.i
    val f: (Int, Int) => Double = (i, j) => {
      if (i == varMinus(vi) && j == varPlus(vi)) {
        val p = toInterval.linearEvaluation(lf)
        2 * math.max(p _1, p _2)
      } else if (i == varPlus(vi) && j == varMinus(vi)) {
        val p = toInterval.linearEvaluation(lf)
        - 2 * math.max(p _1, p _2)
      } else {
        val g1: VarIndex => Boolean = other =>
        vi != other && ((i == varPlus(other) && j == varPlus(vi)) ||
          (i == varMinus(vi) && j == varMinus(other)))
        val g2: VarIndex => Boolean = other =>
        vi != other && ((i == varMinus(other) && j == varPlus(vi)) ||
          (i == varMinus(vi) && j == varPlus(other)))
        val g3: VarIndex => Boolean = other =>
        vi != other && ((i == varPlus(vi) && j == varPlus(other)) ||
          (i == varMinus(other) && j == varMinus(vi)))
        val g4: VarIndex => Boolean = other =>
        vi != other && ((i == varPlus(other) && j == varMinus(vi)) ||
          (i == varPlus(vi) && j == varMinus(other)))

        val chooser = forSomeVar((0 until dimension).map(VarIndex)) _
        val r = (chooser(g1), chooser(g2), chooser(g3), chooser(g4)) match {
          case (Some(other), _, _, _) => {
            val p = toInterval.linearEvaluation(lf - varLf(other, dimension))
            Some(math.max(p _1, p _2))
          }
          case (_, Some(other), _, _) => {
            val p = toInterval.linearEvaluation(lf + varLf(other, dimension))
            Some(math.max(p _1, p _2))
          }
          case (_, _, Some(other), _) => {
            val p = toInterval.linearEvaluation(varLf(other, dimension) - lf)
            Some(math.max(p _1, p _2))
          }
          case (_, _, _, Some(other)) => {
            val p = toInterval.linearEvaluation(- lf - varLf(other, dimension))
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
      case Some(ConstExact(const)) => (m) =>
        e.incrementalClosure(v)(
          singleConstantExactAssignment(v, const.toDouble)(m).elem)
      case Some(SingleExact(Positive, const)) =>
        singlePositiveExactAssignment(v, const.toDouble)
      case Some(SingleExact(Negative, const)) =>
        singleNegativeExactAssignment(v, const.toDouble)
      case Some(DoubleExact(other, Positive, const)) => (matrix) =>
        e.incrementalClosure(v)(
          doublePositiveExactAssignment(v, other, const.toDouble)(matrix).elem)
      case Some(DoubleExact(other, Negative, const)) => (matrix) =>
        e.incrementalClosure(v)(
          doubleNegativeExactAssignment(v, other, const.toDouble)(matrix).elem)
      case None => (matrix) =>
        e.incrementalClosure(v)(thruIntervals(v, lf, dimension)(matrix).elem)
    }
    AbstractOctagon(f(dbm), e)
  }

  sealed trait ExactLinearForm
  case class ConstExact(const: Rational) extends ExactLinearForm
  case class SingleExact(varCoeff: OctaVarCoeff, const: Rational) extends ExactLinearForm
  case class DoubleExact(other: VarIndex, varCoeff: OctaVarCoeff, const: Rational) extends ExactLinearForm

  def nonDeterministicAssignment(n: Int): AbstractOctagon[M] =
    forget(VarIndex(n))

  def linearAssignment(n: Int, lf: LinearForm): AbstractOctagon[M] =
    assignment(VarIndex(n), lf)

  def linearInequality(lf: LinearForm): AbstractOctagon[M] =
    e.decideState(linearInequalityEx(lf).elem) match {
      case CIxed(closed) => AbstractOctagon(closed, e)
      case NCIxed(nc) => AbstractOctagon(e.strongClosure(nc), e)
    }

  // TODO: maybe there's a better option?
  def minimize(lf: LinearForm): RationalExt = toInterval.minimize(lf)

  // TODO: maybe there's a better option?
  def maximize(lf: LinearForm): RationalExt = toInterval.maximize(lf)

  // TODO: maybe there's a better option?
  def frequency(lf: LinearForm): Option[Rational] = toInterval.frequency(lf)

  def cleanup[A](l: List[Option[A]]): List[A] = l match {
    case Nil => Nil
    case Some(x) :: xs => x :: cleanup(xs)
    case None :: xs => cleanup(xs)
  }

  // Every octagonal constraint in a strongly closed DBM defines a half-space
  // that actually touches the octagon. [Mine06]
  def constraints: Seq[LinearForm] = {
    def temp: Double => RationalExt = ???
    val l: List[Option[LinearForm]] =
      (for {
        i <- 0 until e.nOfVars(dbm) * 2
        j <- 0 until e.nOfVars(dbm) * 2
      } yield octaConstrAt(i, j, dbm)
        .map((c) => mapConstraint(temp)(c))
        .flatMap((c) => constrToLf(dimension)(c)))
        .toList

    cleanup(l)
  }

  def isPolyhedral: Boolean = true

  type Domain = OctagonDomain

  private def fromExDBM(eDBM: ExistsMDouble): AbstractOctagon[M] =
    e.decideState(eDBM.elem) match {
      case CIxed(closed) => AbstractOctagon(closed, e)
      case NCIxed(nclosed) => AbstractOctagon(e.strongClosure(nclosed), e)
    }

  def addVariable(): AbstractOctagon[M] = AbstractOctagon(e.addVariable(dbm), e)

  def isTop: Boolean = e.isTopDBM(dbm)
  def isBottom: Boolean = e.isBottomDBM(dbm)

  def delVariable(v: Int): AbstractOctagon[M] =
    AbstractOctagon(e.deleteVariable(VarIndex(v))(dbm), e)

  def mapVariables(rho: Seq[Int]): AbstractOctagon[M] = {
    def converted: VarIndex => Option[VarIndex] = (vi) =>
      if (vi.i >= rho.length || rho(vi.i) == -1) None
      else Some(VarIndex(rho(vi.i)))
    AbstractOctagon(e.mapVariables(converted)(dbm), e)
  }

  def mkString(vars: Seq[String]): String = {
    val sss: IndexedSeq[Option[String]] = for {
      i <- 0 until 2 * e.nOfVars(dbm)
      j <- 0 until 2 * e.nOfVars(dbm)
    } yield octaConstrAt(i, j, dbm).map(prettyConstraint(_.toString, vars))
    cleanup(sss.toList).fold("")((x, y) => x + " ; " + y)
  }

  def isEmpty = isBottom

  def tryCompareTo[B >: AbstractOctagon[M]](that: B)(implicit evidence$1: (B) => PartiallyOrdered[B]): Option[Int] =
    that match {
      case other: AbstractOctagon[M] => {
        e.compare[Double](
          MkEx[Closed, ExistsM](dbm),
          MkEx[Closed, ExistsM](other.dbm)).map {
          case EQ => 0
          case LT => -1
          case GT => 1
        }
      }
    }
}
