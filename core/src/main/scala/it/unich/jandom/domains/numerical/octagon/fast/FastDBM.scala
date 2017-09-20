package it.unich.jandom.domains.numerical.octagon.fast

import breeze.numerics.pow
import it.unich.jandom.domains.numerical.octagon._
import variables._
import CountOps._
import scala.language.higherKinds
import scalaz.{Applicative, Apply, Monoid, Traverse}
import scalaz.std.option._
import scalaz.std.list._
import VarIndexOps._, CountOps._

////////////////////////////////////////////////////////////////////////////////
// FastDBM

sealed trait FastDBM[M[_], SM[_], A] {
  def strongClosure(implicit mev: MEvidence[M, SM], ifield: InfField[A]):
      CFastDBM[M, SM, Closed, A]

  def incrementalClosure(v: VarIndex)(implicit ifield: InfField[A]):
      CFastDBM[M, SM, Closed, A]

  def toFull: FullDBM[M, SM, A]

  def intersection(other: FastDBM[M, SM, A])(implicit ifield: InfField[A]):
      FastDBM[M, SM, A]

  def union(other: FastDBM[M, SM, A])(implicit ifield: InfField[A]):
      FastDBM[M, SM, A]

  def widening(other: FastDBM[M, SM, A])(implicit ifield: InfField[A]):
      FastDBM[M, SM, A]

}

case class MEvidence[M[_], SM[_]](
  dec: Decomposable[M, SM], ds: DenseSparse[M], sub: SubMatrix[SM])

// Full DBMs are fast DBMs that are not decomposed, i.e., they can be either
// dense or sparse.
// Sparsity details, including when to switch between dense and sparse
// representation, is supposed to be handled by the specific implementation of
// the the DenseSparse trait/typeclass.
// An even better thing to do (time permitting) could be to use a suitable
// abstract trait of DBMs that does not talk about sparsity at all (who cares
// if the implementations use a dense/sparse representation anyway, as long as
// they provide the needed DBM-like operations?)
case class FullDBM[M[_], SM[_], A](dbm: M[A], mev: MEvidence[M, SM])
    extends FastDBM[M, SM, A] {

  def toFull = this

  def performStrongClosure(implicit mev: MEvidence[M, SM], ifield: InfField[A]):
      CFastDBM[M, SM, Closed, A] =
    mev.ds.strongClosure(dbm)(ifield) match {
      case Some(closed) => CFast(FullDBM(closed, mev))
      case None         => BottomFast(mev.ds.nOfVars(dbm))
    }

  def strongClosure(implicit mev: MEvidence[M, SM], ifield: InfField[A]):
      CFastDBM[M, SM, Closed, A] = {
    val indepComponents: Seq[Seq[VarIndex]] =
      FastDBMUtils.calculateComponents(dbm, mev.ds)

    if (FastDBMUtils.nuffDecomposed(indepComponents, mev.ds.nOfVars(dbm)))
      DecomposedDBM(dbm, indepComponents, mev).performStrongClosure(mev, ifield)
    else this.performStrongClosure(mev, ifield)
  }

  // Incremental closures are supposed to take advantage of the few modified
  // variables to perform a lightweight closure. We therefore do not try to
  // switch to a decomposed representation.
  def incrementalClosure(v: VarIndex)(implicit ifield: InfField[A])
      : CFastDBM[M, SM, Closed, A] =
    mev.ds.incrementalClosure(v)(dbm) match {
      case Some(m) => CFast(FullDBM(m, mev))
      case None    => BottomFast(mev.ds.nOfVars(dbm))
    }

  def intersection(other: FastDBM[M, SM, A])(implicit ifield: InfField[A]):
      FastDBM[M, SM, A] =
    other match {
      case FullDBM(otherDbm, _) => FullDBM(mev.ds.dbmIntersection(dbm, otherDbm), mev)
      case DecomposedDBM(completeDBM, _, _) =>
        FullDBM(mev.ds.dbmIntersection(dbm, completeDBM), mev)
    }

  def union(other: FastDBM[M, SM, A])(implicit ifield: InfField[A]):
      FastDBM[M, SM, A] = other match {
    case FullDBM(otherDbm, _) => FullDBM(mev.ds.dbmUnion(dbm, otherDbm), mev)
    case DecomposedDBM(completeDBM, _, _) =>
      FullDBM(mev.ds.dbmUnion(dbm, completeDBM), mev)
  }

  def widening(other: FastDBM[M, SM, A])(implicit ifield: InfField[A]):
      FastDBM[M, SM, A] = other match {
    case FullDBM(otherDbm, _) => FullDBM(mev.ds.widening(dbm, otherDbm), mev)
    case DecomposedDBM(completeDBM, _, _) =>
      FullDBM(mev.ds.widening(dbm, completeDBM), mev)
  }

}

// We store the independent components as a linked list of linked lists of
// variable indices.
// Independent components correspond to submatrices in the complete DBM,
// and these can be dense or sparse. The exact sparsity is computed
// on-the-fly before closure.

// Octagon operators on the decomposed type:
// Apart from closure (which is specialized according to sparsity), we apply the
// standard operators independently on each submatrix.

// NNI for decomposed DBMs is computed as the sum of its submatrices' NNI.
case class DecomposedDBM[M[_], SM[_], A](completeDBM: M[A],
                                         indepComponents: Seq[Seq[VarIndex]],
                                         mev: MEvidence[M, SM]) extends FastDBM[M, SM, A] {

  def toFull: FullDBM[M, SM, A] = FullDBM(completeDBM, mev)

  def performStrongClosure(implicit mev: MEvidence[M, SM], ifield: InfField[A]):
      CFastDBM[M, SM, Closed, A] = {

    val submatrices = indepComponents.map(seq => mev.dec.extract(seq)(completeDBM))
    val closed = submatrices.map(m => mev.sub.strongClosure(m))

    Applicative[Option].sequence(closed.toList) match {

      case Some(closedSubs) => {
        val newMatrix = closedSubs.foldRight(completeDBM)(
          (sub, full) => mev.dec.pour(sub)(full))
        CFast(DecomposedDBM(newMatrix, indepComponents, mev))
      }

      case None => BottomFast(mev.ds.nOfVars(completeDBM))
    }

  }

  def strongClosure(implicit mev: MEvidence[M, SM], ifield: InfField[A]):
      CFastDBM[M, SM, Closed, A] = {
    val indepComponents: Seq[Seq[VarIndex]] =
      FastDBMUtils.calculateComponents(completeDBM, mev.ds)

    if (FastDBMUtils.nuffDecomposed(indepComponents, mev.ds.nOfVars(completeDBM)))
      DecomposedDBM(completeDBM, indepComponents, mev).performStrongClosure(mev, ifield)
    else toFull.performStrongClosure(mev, ifield)
  }

  def incrementalClosure(v: VarIndex)(implicit ifield: InfField[A]):
      CFastDBM[M, SM, Closed, A] =

    indepComponents.find(_.contains(v)) match {
      case Some(comp) => {
        val subMat = mev.dec.extract(comp)(completeDBM)
        mev.sub.incrementalClosure(v)(subMat) match {
          case Some(closed) => {
            val newMat = mev.dec.pour(closed)(completeDBM)
            CFast(DecomposedDBM(newMat, indepComponents, mev))
          }
          case None => BottomFast(mev.ds.nOfVars(completeDBM))
        }
      }
      case None => CFast(this)
    }

  def intersection(other: FastDBM[M, SM, A])(implicit ifield: InfField[A]):
      FastDBM[M, SM, A] =
    other match {

      case FullDBM(otherFull, _) =>
        FullDBM(mev.ds.dbmIntersection(completeDBM, otherFull), mev)

      case DecomposedDBM(otherCompleteDBM, otherComponents, _) =>

        // TODO: maybe too approximate?

        // Compute the union of all independent components, yielding the set of
        // all variables that can be different from infinity (assuming that
        // all operators preserve soundness of the components index.)
        val vars1 = indepComponents.foldRight(Seq[VarIndex]())(_ ++ _).toSet
        val vars2 = otherComponents.foldRight(Seq[VarIndex]())(_ ++ _).toSet
        val component = (vars1 union vars2).toSeq

        // create new submatrices with the same components
        val sub1 = mev.dec.extract(component)(completeDBM)
        val sub2 = mev.dec.extract(component)(otherCompleteDBM)

        val newMat = mev.sub.dbmIntersection(sub1, sub2)
        DecomposedDBM(mev.dec.pour(newMat)(completeDBM), Seq(component), mev)
    }

  def union(other: FastDBM[M, SM, A])(implicit ifield: InfField[A]):
      FastDBM[M, SM, A] = other match {

    case FullDBM(otherFull, _) =>
      FullDBM(mev.ds.dbmUnion(completeDBM, otherFull), mev)

    case DecomposedDBM(otherCompleteDBM, otherComponents, _) =>

      val vars1 = indepComponents.foldRight(Seq[VarIndex]())(_ ++ _).toSet
      val vars2 = otherComponents.foldRight(Seq[VarIndex]())(_ ++ _).toSet
      val vars = vars1 intersect vars2
      val newComps = indepComponents.map(_.filter(vars.contains(_)))
      val matrices = newComps.map(c => {
        val sub1 = mev.dec.extract(c)(completeDBM)
        val sub2 = mev.dec.extract(c)(otherCompleteDBM)
        mev.sub.dbmUnion(sub1, sub2) })
      val newMat = matrices.foldLeft(completeDBM)(
        (mat, subMat) => mev.dec.pour(subMat)(mat))
      DecomposedDBM(newMat, newComps, mev)
  }

  def widening(other: FastDBM[M, SM, A])(implicit ifield: InfField[A]):
      FastDBM[M, SM, A] = other match {

    case FullDBM(otherFull, _) =>
      FullDBM(mev.ds.widening(completeDBM, otherFull), mev)

    case DecomposedDBM(otherCompleteDBM, otherComponents, _) =>
      val vars1 = indepComponents.foldRight(Seq[VarIndex]())(_ ++ _).toSet
      val vars2 = otherComponents.foldRight(Seq[VarIndex]())(_ ++ _).toSet
      val vars = vars1 intersect vars2
      val newComps = indepComponents.map(_.filter(vars.contains(_)))
      val matrices = newComps.map(c => {
        val sub1 = mev.dec.extract(c)(completeDBM)
        val sub2 = mev.dec.extract(c)(otherCompleteDBM)
        mev.sub.widening(sub1, sub2) })
      val newMat = matrices.foldLeft(completeDBM)(
        (mat, subMat) => mev.dec.pour(subMat)(mat))
      DecomposedDBM(newMat, newComps, mev)
  }

}

////////////////////////////////////////////////////////////////////////////////
// Closable FastDBMs

// ADT of "closable" DBMs in their fast implementation from Vechev et al.
// They are "closable" in the sense that they augment the ADT of fast DBMs with
// the type-level capability of being indexed by their strong closure state.
sealed trait CFastDBM[M[_], SM[_], _, A]

// Constructor of *closed* fast DBMs.
case class CFast[M[_], SM[_], A](m: FastDBM[M, SM, A])
    extends CFastDBM[M, SM, Closed, A]

// Constructor of *not-necessarily-closed* fast DBMs.
case class NCFast[M[_], SM[_], A](m: FastDBM[M, SM, A])
    extends CFastDBM[M, SM, NonClosed, A]

case class BottomFast[M[_], SM[_], A](nOfVars: VarCount)
    extends CFastDBM[M, SM, Closed, A]

////////////////////////////////////////////////////////////////////////////////

object FastDBMUtils {

  // TODO: find a better definition.
  // For example, consider nuff decomposed those components
  // with size <= 50% of the total.
  def nuffDecomposed(indepComps: Seq[Seq[VarIndex]], count: VarCount): Boolean =
    indepComps.length >= 2

  def fastInnerMatrix[M[_], SM[_], S <: DBMState, A](fdbm: FastDBM[M, SM, A]): M[A] =
    fdbm match {
      case FullDBM(m: M[A], _) => m
      case DecomposedDBM(m: M[A], _, _) => m
    }

  def nOfVars[M[_], SM[_], S, A](dbm: CFastDBM[M, SM, S, A])
                         (implicit mev: MEvidence[M, SM]): VarCount =
    dbm match {
      case CFast(m: FastDBM[M, SM, A]) => mev.ds.nOfVars(fastInnerMatrix(m))
      case NCFast(m: FastDBM[M, SM, A]) => mev.ds.nOfVars(fastInnerMatrix(m))
      case BottomFast(n) => n
    }

  def calculateComponents[M[_], A](m: M[A], ds: DenseSparse[M])
    (implicit ifield: InfField[A]): List[List[VarIndex]] = {

    def related(vi: VarIndex, vj: VarIndex): Boolean =
      Set((varPlus(vi), varPlus(vj)), (varPlus(vi), varMinus(vj)),
          (varMinus(vi), varPlus(vj)), (varMinus(vi), varMinus(vj)))
        .filter({ case (i, j) => i != j})
        .exists({ case (i, j) => ifield.!=(ds.get(i, j)(m), ifield.infinity)})

    val nOfVars = ds.nOfVars(m)
    val rels =
      for (vi <- allVars(nOfVars); vj <- allVars(nOfVars); if related(vi, vj))
      yield (vi, vj)

    val indComps =
      rels.foldLeft(Set[Set[VarIndex]]())({ case (comps, (vi, vj)) =>
        val compI = comps.find(_.contains(vi)).getOrElse(Set())
        val compJ = comps.find(_.contains(vj)).getOrElse(Set())
        val newComp = Set(vi, vj) ++ compI ++ compJ
        comps.filter(c => !c.contains(vi) && !c.contains(vj)) + newComp
      })

    indComps.map(_.toList).toList
  }

}
