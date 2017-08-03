package it.unich.jandom.domains.numerical.octagon

import it.unich.jandom.domains.numerical._

import scala.language.higherKinds
import scala.language.postfixOps
import VarIndexOps._
import VarIndexUtils._
import spire.math.Rational
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

case class AbstractOctagon[D <: NumericalDomain, M[_, _]](
  dbm: M[Closed, Double],
  d: D, e: DifferenceBoundMatrix[M] { type PosetConstraint[A] = InfField[A] })
    extends NumericalProperty[AbstractOctagon[D, M]] {

  import AbstractOctagon._
  import DBMUtils._

  private def withDBM(dbm: M[Closed, Double]): AbstractOctagon[D, M] =
    AbstractOctagon(dbm, d, e)

  def dimension: Int = e.nOfVars(dbm)

  def union(other: AbstractOctagon[D, M]): AbstractOctagon[D, M] =
    withDBM(e.dbmUnion(dbm, other.dbm)(InfField.infFieldDouble))

  def intersection(other: AbstractOctagon[D, M]): AbstractOctagon[D, M] =
    withDBM(e.strongClosure(e.dbmIntersection(dbm, other.dbm).elem))

  def forget(vi: VarIndex): AbstractOctagon[D, M] = withDBM(e.forget(vi)(dbm))

  def top = withDBM(e.topDBM[Double](e.nOfVars(dbm)))
  def bottom = withDBM(e.bottomDBM[Double](e.nOfVars(dbm)))

  def widening(other: AbstractOctagon[D, M]): AbstractOctagon[D, M] =
    withDBM(e.strongClosure(e.widening(dbm, other.dbm).elem))

  def narrowing(other: AbstractOctagon[D, M]): AbstractOctagon[D, M] =
    withDBM(e.strongClosure(e.narrowing(dbm, other.dbm).elem))

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
    val isEmpty : Boolean =
      low.forall(_ == Double.PositiveInfinity) &&
      high.forall(_ == Double.NegativeInfinity)
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
  private [numerical] def fallbackUpdate (lf: LinearForm) : ExistsDBM[({ type T[S] = M[S, Double]})#T] = {
    def inverseVarPlus = inverseVarPlusMinus(VarIndexOps.Positive, this.dimension * 2, _: Int)
    def inverseVarMinus = inverseVarPlusMinus(VarIndexOps.Negative, this.dimension * 2, _: Int)
    /**
      * Computes enhanced fallback abstraction according to Mine06 fig. 21, ie:
      *
      * ({ e<= 0?}(m))_{ij} = min(m_{ij}, m'_{ij})
      *
      * The cases for m'_{ij} are:
      *
      *  1.  2*max (Int (V_j0 - e))
      *  2. -2*max (int(-V_j0 - e))
      *  3a. max(Int(V_j0 - V_i0 - e))
      *  3b.    "   "   "
      *  4.  max(Int( V_j0 + V_i0- e))
      *  5.  max(Int(-V_j0 - V_i0- e))
      *  6. id
      *
      * Visually the cases for fallback map as follows:
      *
      *     |+ |- |+ |- |+ |- |
      *   --+=====+--+--+--+--+
      *    +|id|1 |3b|4 |3b|4 |
      *   --+--+--+--+--+--+--+
      *    -| 2|id|5 |3a|5 |3a|
      *   --+=====+=====+--+--+
      *    +|3b|4 |id|1 |3b|4 |
      *   --+--+--+--+--+--+--+
      *    -|5 |3a| 2|id|5 |3a|
      *   --+--+--+=====+=====+
      *    +|3b|4 |3b|1 |id|1 |
      *   --+--+--+--+--+--+--+
      *    -|5 |3a|5 |3a| 2|id|
      *   --+--+--+--+--+=====+
      */
    def g (i: Int, j: Int) : Double = {
      if (inverseVarPlus(i) != None & inverseVarMinus(j) != None
        & inverseVarPlus(i) == inverseVarMinus(j)
      ) { // Case 1 with i/2 = j0
        val j0 = inverseVarPlus(i).get
        2*(this.toInterval.maximize(LinearForm.v(j0.i) - lf)).toDouble
     } else if (inverseVarMinus(i) != None & inverseVarPlus(j) != None
       & inverseVarMinus(i) == inverseVarPlus(j)
      ) { // Case 2 with j/2 = j0
        val j0 = inverseVarPlus(j).get
        -2*(this.toInterval.maximize(LinearForm.v(j0.i) - lf)).toDouble
      } else if (inverseVarMinus(j) != None & inverseVarMinus(i) != None
        & inverseVarMinus(i) != inverseVarMinus(j)
      ) { // Case 3a, (i+1)/2 = i0, (j+1)/2 = j0
        val j0 = inverseVarMinus(j).get
        val i0 = inverseVarMinus(i).get
        this.toInterval.maximize(LinearForm.v(j0.i) - LinearForm.v(i0.i) - lf).toDouble
      } else if (inverseVarPlus(j) != None & inverseVarPlus(i) != None
        & inverseVarPlus(i) != inverseVarPlus(j)
      ) { // Case 3b, i/2 = i0, j/2 = j0
        val j0 = inverseVarPlus(j).get
        val i0 = inverseVarPlus(i).get
        this.toInterval.maximize(LinearForm.v(j0.i) - LinearForm.v(i0.i) - lf).toDouble
      } else if (inverseVarMinus(j) != None & inverseVarPlus(i) != None
        & inverseVarMinus(j) != inverseVarPlus(i)
      ) { // Case 4, i0 = i/2, j0 = j+1/2
        val j0 = inverseVarMinus(j).get
        val i0 = inverseVarPlus(i).get
        this.toInterval.maximize(LinearForm.v(j0.i) + LinearForm.v(i0.i) - lf).toDouble
      } else if (inverseVarMinus(i) != None & inverseVarPlus(j) != None
        & inverseVarPlus(j) != inverseVarMinus(i)
      ) { // Case 5, j0 = j, i0 = (i+1)/2
        val j0 = inverseVarPlus(j).get
        val i0 = inverseVarMinus(i).get
        this.toInterval.maximize(- LinearForm.v(j0.i) - LinearForm.v(i0.i) - lf).toDouble
      } else {
        e.get(i,j)(dbm).get
      }
    }
    val f = (i : Int, j : Int) => math.min(
      e.get(i,j)(dbm).get,
          g(i,j)
    )
    e.update(f)(dbm)
  }

  def linearInequalityEx (lf: LinearForm): ExistsDBM[({ type T[S] = M[S, Double]})#T] = {
    val t : AbstractTest = decodeTest(lf)
    t match {
      case Fallback() => fallbackUpdate(lf)
      case et : ExactTest => {
        val f = et match {
          case Case1Test(vj, c) => {
            (i : Int, j : Int) =>
            if (i == varPlus(VarIndex(vj)) & j == varMinus(VarIndex(vj)))
              math.min(-2*c.toDouble, e.get(i,j)(dbm).get)
            else
              e.get(i,j)(dbm).get
          }
          case Case2Test(vj, c) => {
            (i : Int, j : Int) =>
            if (i == varPlus(VarIndex(vj)) & j == varMinus(VarIndex(vj)))
              math.min(-2*c.toDouble, e.get(i,j)(dbm).get)
            else
              e.get(i,j)(dbm).get
          }
          case Case3Test(vj, vi, c) => {
            (i : Int, j : Int) =>
            if ((i == varPlus(VarIndex(vi)) & j == varPlus(VarIndex(vj)))
              | (i == varMinus(VarIndex(vj)) & j == varMinus(VarIndex(vi)))
            )
              math.min(-1*c.toDouble, e.get(i,j)(dbm).get)
            else
              e.get(i,j)(dbm).get
          }
          case Case4Test(vj, vi, c) => {
            (i : Int, j : Int) =>
            if ((i == varMinus(VarIndex(vi)) & j == varPlus(VarIndex(vj)))
              | (i == varMinus(VarIndex(vj)) & j == varPlus(VarIndex(vi)))
            )
              math.min(-1*c.toDouble, e.get(i,j)(dbm).get)
            else
              e.get(i,j)(dbm).get
          }
          case Case5Test(vj, vi, c) => {
            (i : Int, j : Int) =>
            if ((i == varPlus(VarIndex(vi)) & j == varMinus(VarIndex(vj)))
              | (i == varPlus(VarIndex(vj)) & j == varMinus(VarIndex(vi)))
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

  def assignment(v: VarIndex, lf: LinearForm): AbstractOctagon[D, M] = {
    val f: M[Closed, Double] => M[Closed, Double] = decideLinearForm(v, lf) match {
      case Some(ConstExact(const)) => (m) =>
        e.incrementalClosure(v)(
          singleConstantExactAssignment(v, const.toDouble)(m, e).elem)
      case Some(SingleExact(Positive, const)) => (m) =>
        singlePositiveExactAssignment(v, const.toDouble)(m, e)
      case Some(SingleExact(Negative, const)) => (m) =>
        singleNegativeExactAssignment(v, const.toDouble)(m, e)
      case Some(DoubleExact(other, Positive, const)) => (matrix) =>
        e.incrementalClosure(v)(
          doublePositiveExactAssignment(v, other, const.toDouble)(matrix, e).elem)
      case Some(DoubleExact(other, Negative, const)) => (matrix) =>
        e.incrementalClosure(v)(
          doubleNegativeExactAssignment(v, other, const.toDouble)(matrix, e).elem)
      case None => (matrix) =>
        e.incrementalClosure(v)(thruIntervals(v, lf, dimension)(matrix).elem)
    }
    withDBM(f(dbm))
  }

  sealed trait OctaVarCoeff
  object Positive extends OctaVarCoeff { }
  object Negative extends OctaVarCoeff { }
  sealed trait ExactLinearForm
  case class ConstExact(const: Rational) extends ExactLinearForm
  case class SingleExact(varCoeff: OctaVarCoeff, const: Rational) extends ExactLinearForm
  case class DoubleExact(other: VarIndex, varCoeff: OctaVarCoeff, const: Rational) extends ExactLinearForm

  def nonDeterministicAssignment(n: Int): AbstractOctagon[D, M] =
    forget(VarIndex(n))

  def linearAssignment(n: Int, lf: LinearForm): AbstractOctagon[D, M] =
    assignment(VarIndex(n), lf)

  def linearInequality(lf: LinearForm): AbstractOctagon[D, M] =
    e.decideState(linearInequalityEx(lf).elem) match {
      case CIxed(closed) => withDBM(closed)
      case NCIxed(nc) => withDBM(e.strongClosure(nc))
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
      } yield octaConstrAt(i, j, dbm)(InfField.infFieldDouble, e)
        .map((c) => mapConstraint(temp)(c))
        .flatMap((c) => constrToLf(dimension)(c)))
        .toList

    cleanup(l)
  }

  def isPolyhedral: Boolean = true

  type Domain = D

  private def fromExDBM(eDBM: ExistsDBM[({ type T[S] = M[S, Double]})#T]): AbstractOctagon[D, M] =
    e.decideState(eDBM.elem) match {
      case CIxed(closed) => withDBM(closed)
      case NCIxed(nclosed) => withDBM(e.strongClosure(nclosed))
    }

  def addVariable(): AbstractOctagon[D, M] = withDBM(e.addVariable(dbm))

  def isTop: Boolean = e.isTopDBM(dbm)
  def isBottom: Boolean = e.isBottomDBM(dbm)

  def delVariable(v: Int): AbstractOctagon[D, M] =
    withDBM(e.deleteVariable(VarIndex(v))(dbm))

  def mapVariables(rho: Seq[Int]): AbstractOctagon[D, M] = {
    def converted: VarIndex => Option[VarIndex] = (vi) =>
      if (vi.i >= rho.length || rho(vi.i) == -1) None
      else Some(VarIndex(rho(vi.i)))
    withDBM(e.mapVariables(converted)(dbm))
  }

  def mkString(vars: Seq[String]): String = {
    val sss: IndexedSeq[Option[String]] = for {
      i <- 0 until 2 * e.nOfVars(dbm)
      j <- 0 until 2 * e.nOfVars(dbm)
    } yield octaConstrAt(i, j, dbm)(InfField.infFieldDouble, e)
      .map(prettyConstraint(_.toString, vars))
    cleanup(sss.toList).fold("")((x, y) => x + " ; " + y)
  }

  def domain: this.Domain = d

  def isEmpty = isBottom

  type ExistsM[S] = M[S, Double]

  def tryCompareTo[B >: AbstractOctagon[D, M]](that: B)(implicit evidence$1: (B) => PartiallyOrdered[B]): Option[Int] =
    that match {
      case other: AbstractOctagon[D, M] => {
        e.compare[Double](
          MkEx[Closed, ExistsM](dbm),
          MkEx[Closed, ExistsM](other.dbm)).map {
          case EQ => 0
          case LT => -1
          case GT => 1
        }
      }
    }
  def thruIntervals(vi: VarIndex, lf: LinearForm, dimension: Int)
    (dbm: M[Closed, Double]): ExistsDBM[({ type T[S] = M[S, Double]})#T] = {
    // val v = vi.i
    val f: (Int, Int) => Double = (i, j) => {
      if (i == varMinus(vi) && j == varPlus(vi)) {
        val p = toInterval.linearEvaluation(lf)
        2 * math.max(p _1, p _2)
      } else if (i == varPlus(vi) && j == varMinus(vi)) {
        val p = toInterval.linearEvaluation(lf)
        - 2 * math.min(p _1, p _2)
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
          case None => (e.get(i,j)(e.strongClosure(dbm))).get
        }
      }
    }
    e.update(f)(dbm)
  }

  private def varLf(v: VarIndex, dimension: Int): LinearForm = {
    val sss: Seq[Rational] =
      (0 to dimension).map(x => if (x == v.i + 1) Rational(1) else Rational(0))
    new DenseLinearForm(sss)
  }

  // Evaluation of linear assignment using interval arithmetics.
  def lfAsInterval(v: VarIndex, lf: LinearForm): (Double, Double) = ???
}

sealed abstract class AbstractTest
case class Fallback() extends AbstractTest
sealed abstract class ExactTest extends AbstractTest
case class Case1Test(val vj : Int, val c : Rational) extends ExactTest
case class Case2Test(val vj : Int, val c : Rational) extends ExactTest
case class Case3Test(val vj : Int, val vi : Int, val c : Rational) extends ExactTest
case class Case4Test(val vj : Int, val vi : Int, val c : Rational) extends ExactTest
case class Case5Test(val vj : Int, val vi : Int, val c : Rational) extends ExactTest

object AbstractOctagon {
  /**
    * Decodes a Linearform into an AbstractTest, ie one of 6 cases discussed in Mine06.
    */
  private[octagon] def decodeTest (lf : LinearForm) : AbstractTest = {
    require (lf.coeffs.size > 0)
    if ( lf.homcoeffs.exists { !List(1,0,-1).contains(_) }
       | lf.homcoeffs.filter(_ != 0).size > 2
       | lf.homcoeffs.filter(_ != 0).size == 0
       // TODO: perhaps case (c,0,0...,0) needs is a special case of one of the others?
    ) {
      Fallback()
      // Mine06 fig 20 does not give an exact abstraction for > 2
      // coeffs != 0 or any coeffs != 1, -1
    } else {
      val c = lf.known
      val nonZeroPairs = lf.pairs.filter(_._2 != 0)
      if (nonZeroPairs.size == 1) {
        // Cases 1, 2
        val vl = nonZeroPairs.head
        if (vl._2 == 1)
          Case1Test (vl._1, c)
        else {
          require(vl._2 == -1)
          Case2Test (vl._1, c)
        }
      } else if (nonZeroPairs.size == 2) {
        if (nonZeroPairs.exists(_._2 == 1)) {
          // Cases 3, 4
          val vl = nonZeroPairs.filter(_._2 == 1).head
          val vk = nonZeroPairs.diff(List(vl)).head
          if (vk._2 == -1)
            Case3Test (vl._1, vk._1, c)
          else
            Case4Test (vl._1, vk._1, c)
        } else {
          // Case 5
          val vl = nonZeroPairs.head
          val vk = nonZeroPairs.diff(List(vl)).head
          Case5Test (vl._1, vk._1, c)
        }
      } else {
        // TODO: temporary, eventually remove
        throw new RuntimeException("If you got here my logic is broken, sorry...")
      }
    }
  }

  def fromInterval[D <: NumericalDomain, M[_,_]]
    (box: BoxDoubleDomain#Property,
      d: D, e: DifferenceBoundMatrix[M] { type PosetConstraint[A] = InfField[A] })
      : AbstractOctagon[D, M] = {
    require (box.low.size == box.high.size)
    val indices = (0 until box.high.size).map(x => VarIndex(x))
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
    AbstractOctagon(e.fromFun(box.dimension*2, f), d, e)
  }
}

object DBMUtils {
  type ExistsMDouble[M[_,_]] = ExistsDBM[({ type T[S] = M[S, Double]})#T]

  def singleConstantExactAssignment[M[_,_], S <: DBMState]
    (v: VarIndex, const: Double)
    (dbm: M[S, Double], e: DifferenceBoundMatrix[M]): ExistsMDouble[M] = {
    val f: (Int, Int) => Double = (i, j) =>
      if (i == varPlus(v) && j == varMinus(v)) -2 * const else
        if (i == varMinus(v) && j == varPlus(v)) 2 * const else
          (e.get(i, j)(e.forget(v)(dbm))).get
    e.update(f)(dbm)
  }

  // single exact assignments preserve strong closure
  def singlePositiveExactAssignment[M[_,_], S <: DBMState]
    (v: VarIndex, const: Double)
    (dbm: M[S, Double], e: DifferenceBoundMatrix[M]) : M[S, Double] = e.addScalarOnVar(v, const)(dbm)

  def doublePositiveExactAssignment[M[_,_], S <: DBMState]
    (vi: VarIndex, vother: VarIndex, const: Double)
    (dbm: M[S, Double], e: DifferenceBoundMatrix[M]): ExistsMDouble[M] = {

    val f: (Int, Int) => Double = (i, j) => {
      val g1 = (i == varPlus(vi) && j == varPlus(vother)) ||
               (i == varMinus(vother) && j == varMinus(vi))
      val g2 = (i == varPlus(vother) && j == varPlus(vi)) ||
               (i == varMinus(vi) && j == varMinus(vother))
      if (g1) -const else
        if (g2) const else
          (e.get(i, j)(e.forget(vi)(e.strongClosure(dbm)))).get
    }
    e.update(f)(dbm)
  }

  // x := - x
  // this preserves strong closure
  def singleNegativeZeroAssignment[M[_,_], S <: DBMState]
    (v: VarIndex)(dbm: M[S, Double], e: DifferenceBoundMatrix[M]): M[S, Double] = e.flipVar(v)(dbm)

  // x := - y
  def doubleNegativeZeroAssignment[M[_,_], S <: DBMState]
    (v: VarIndex, other: VarIndex)(dbm: M[S, Double], e: DifferenceBoundMatrix[M]): ExistsMDouble[M] = {
    val m = doublePositiveExactAssignment(v, other, 0)(dbm, e)
    val mm = singleNegativeZeroAssignment(v)(m.elem, e)
    MkEx[m.State, ({ type T[A] = M[A, Double] })#T](mm)
  }

  // x := - x + c
  def singleNegativeExactAssignment[M[_,_], S <: DBMState] (v: VarIndex, const: Double)
                                   (dbm: M[Closed, Double],  e: DifferenceBoundMatrix[M]): M[Closed, Double] =
    singlePositiveExactAssignment(v, const)(
      singleNegativeZeroAssignment(v)(dbm, e), e)

  // x := - y + c
  def doubleNegativeExactAssignment[M[_,_], S <: DBMState]
    (v: VarIndex, other: VarIndex, const: Double)
    (dbm: M[S, Double],  e: DifferenceBoundMatrix[M]): ExistsMDouble[M] = {
    val m = doubleNegativeZeroAssignment(v, other)(dbm, e)
    val mm = singlePositiveExactAssignment(v, const)(m.elem, e)
    MkEx[m.State, ({ type T[A] = M[A, Double] })#T](mm)
  }
}
