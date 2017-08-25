package it.unich.jandom.domains.numerical.octagon

import it.unich.jandom.domains.numerical._

import scala.language.higherKinds
import scala.language.postfixOps
import VarIndexOps._
import VarIndexUtils._
import it.unich.jandom.domains.numerical.octagon.CountOps.allIndices
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

  def dimension: Int = e.nOfVars(dbm).count

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

  import CountOps._

  def toInterval: BoxDoubleDomain#Property = {
    val closed = e.strongClosure(dbm)
    val l: List[(Double, Double)] =
      allVars(e.nOfVars(dbm)).map(v => projectInterval(v, closed)).toList
    val (low, high) = l.unzip
    createInterval(low.toArray, high.toArray, isEmpty = false)
  }

  def linearInequalityEx (lf: LinearForm): ExistsDBM[({ type T[S] = M[S, Double]})#T] = ???
  //////////////////////// BOOL TESTS //////////////////////////////////////////

  def assignment(v: VarIndex, lf: LinearForm): AbstractOctagon[D, M] = ???
  /**
   * Computes intersection with `lf != 0`.
   *
   * Amounts to identity as Mine06 does not give a better abstraction
   * than id
   */
  def linearDisequality (lf: LinearForm) = this

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
        i <- allIndices(doubledVarCount(e.nOfVars(dbm)))
        j <- allIndices(doubledVarCount(e.nOfVars(dbm)))
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
    val indexedVars = vars.zipWithIndex
                          .map({ case (name, i) => (name, VarIndex(i))})
    val str =
      for ((varj, vj) <- indexedVars;
           (namej, posj) <- Seq((s"+$varj", varPlus(vj)),
                                (s"-$varj", varMinus(vj)));
           (vari, vi) <- indexedVars;
           (namei, posi) <- Seq((s"+$vari", varPlus(vi)),
                                (s"-$vari", varMinus(vi)))) yield {
        val elem = e.get(posi, posj)(dbm)
        s"$namej - $namei <= $elem"
      }
    "[ " + str.reduceOption(_ ++ " ; " ++ _).getOrElse("") + " ]"
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

object DBMUtils {
  type ExistsMDouble[M[_,_]] = ExistsDBM[({ type T[S] = M[S, Double]})#T]

  def singleConstantExactAssignment[M[_,_], S <: DBMState]
    (v: VarIndex, const: Double)
    (dbm: M[S, Double], e: DifferenceBoundMatrix[M]): ExistsMDouble[M] = {
    val f: (Int, Int) => Double = (i, j) =>
      if (i == varPlus(v) && j == varMinus(v)) 2 * const else
        if (i == varMinus(v) && j == varPlus(v)) -2 * const else
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
