package it.unich.jandom.domains.numerical.octagon.fast

import breeze.numerics.pow
import it.unich.jandom.domains.numerical.octagon._
import variables._
import CountOps._
import scala.language.higherKinds
import scalaz.{Applicative, Apply, Monoid, Traverse}
import scalaz.std.option._
import scalaz.std.list._

// Sparsity D = 1 - (nni / (2n^2 + 2n))
// Switching between DBMs: the sparsity can increase, for instance during
// widening. Recovering sparsity information and  independent components has a
// quadratic worst case complexity, so we only perform it by piggybacking on
// the closure operator. We also use closure computations as switching points.


case class MEvidence[M[_], SM[_]](dec: Decomposable[M, SM],
                                  ds: DenseSparse[M],
                                  sub: SubMatrix[SM])

// DifferenceBoundMatrix instance for CFastDBM.
object CFDBMInstance {

  def instance[M[_], SM[_]](implicit mev: MEvidence[M, SM]) =
    new DifferenceBoundMatrix[({ type T[S, A] = CFastDBM[M, SM, S, A] })#T] {

      type PosetConstraint[A] = InfField[A]

      def update[S <: DBMState, A](f: (Int, Int) => A)(m: CFastDBM[M, SM, S, A])
                                  (implicit ifield: InfField[A]): ExistsM[A] =
        Utils.liftFromInner(mev.ds.update(f))(m)

      def incrementalClosure[S <: DBMState, A](v: VarIndex)
                               (dbm: CFastDBM[M, SM, S, A])
        (implicit evidence: InfField[A]): CFastDBM[M, SM, Closed, A] =
        dbm match {
          case BottomFast(n) => BottomFast(n)
          case CFast(m: FastDBM[M, SM, A]) => CFast(m)
          case NCFast(m: FastDBM[M, SM, A]) => m.incrementalClosure(v)
        }

      def strongClosure[S <: DBMState, A](dbm: CFastDBM[M, SM, S, A])
                          (implicit evidence: InfField[A]): CFastDBM[M, SM, Closed, A] =
        dbm match {
          case BottomFast(n) => BottomFast(n)
          case CFast(m: FastDBM[M, SM, A]) => CFast(m)
          case NCFast(m: FastDBM[M, SM, A]) => m.strongClosure
        }

      def forget[S <: DBMState, A](vi: VarIndex)(m: CFastDBM[M, SM, S, A])
                                  (implicit f: InfField[A]): CFastDBM[M, SM, S, A] =
        Utils.mapFastDBM[M, SM, S, A](Utils.mapInnerMatrix(mev.ds.forget(vi)))(m)

      def nOfVars[S <: DBMState, A](m: CFastDBM[M, SM, S, A]): VarCount = Utils.nOfVars(m)

      def get[S <: DBMState, A](i: Int, j: Int)(m: CFastDBM[M, SM, S, A])
                               (implicit ifield: InfField[A]): Option[A] =
        Utils.cfastInnerMatrix(m).map((inner) => mev.ds.get(i, j)(inner))

      def dbmIntersection[A, R <: DBMState, S <: DBMState]
        (m1: CFastDBM[M, SM, R, A], m2: CFastDBM[M, SM, S, A])
        (implicit ifield: InfField[A]): ExistsM[A] = {

        def aux(m1: FastDBM[M, SM, A], m2: FastDBM[M, SM, A]): FastDBM[M, SM, A] =
          (m1, m2) match {
            case (FullDBM(dbm1, mev), FullDBM(dbm2, _)) =>
              FullDBM(mev.ds.dbmUnion(dbm1, dbm2), mev)
            case (DecomposedDBM(dbm1, comps1, mev), DecomposedDBM(dbm2, comps2, _)) =>
              // compute sets of initialized variables, i.e., variables
              // that can be different from infinity
              val vars1 = comps1.foldRight(Seq[VarIndex]())(_ ++ _).toSet
              val vars2 = comps2.foldRight(Seq[VarIndex]())(_ ++ _).toSet

              // create a new component with the variables that can be
              // different from infinity in the intersection
              val vars = vars1 union vars2
              val component = vars.toSeq

              // create new submatrices with the same components
              val sub1 = mev.dec.extract(component)(dbm1)
              val sub2 = mev.dec.extract(component)(dbm2)

              val newMat = mev.sub.dbmIntersection(sub1, sub2)
              ??? // DecomposedDBM(newMat, Seq(component), mev)
            case (dbm1 @ DecomposedDBM(_, _, _), dbm2 @ FullDBM(_, _)) =>
              aux(dbm1.toFull, dbm2)
            case (dbm1 @ FullDBM(_, _), dbm2 @ DecomposedDBM(_, _, _)) =>
              aux(dbm1, dbm2.toFull)
          }

        (m1, m2) match {
          case (BottomFast(nOfVars), _) =>
            Utils.packEx(BottomFast(nOfVars))
          case (_, BottomFast(nOfVars)) =>
            Utils.packEx(BottomFast(nOfVars))
          case (CFast(dbm1), CFast(dbm2)) =>
            Utils.packEx(NCFast(aux(dbm1, dbm2)))
          case (CFast(dbm1), NCFast(dbm2)) =>
            Utils.packEx(NCFast(aux(dbm1, dbm2)))
          case (NCFast(dbm1), CFast(dbm2)) =>
            Utils.packEx(NCFast(aux(dbm1, dbm2)))
          case (NCFast(dbm1), NCFast(dbm2)) =>
            Utils.packEx(NCFast(aux(dbm1, dbm2)))
        }
      }

      def topDBM[A](nOfVars: VarCount)(implicit ifield: InfField[A]):
          CFastDBM[M, SM, Closed, A] = {
        val top = mev.dec.pure(nOfVars, ifield.infinity)
        val cTop = allIndices(varCountToDim(nOfVars))
          .foldLeft(top)((m, i) => mev.ds.update(i, i, ifield.zero)(m))
        CFast(FullDBM(cTop, mev))
      }

      def bottomDBM[A](nOfVars: VarCount)(implicit ifield: InfField[A]): CFastDBM[M, SM, Closed, A] =
        BottomFast(nOfVars)

      def fromFun[A](d: Dimension, f: ((Int, Int) => A))(implicit ifield: InfField[A]): CFastDBM[M, SM, Closed, A] =
        CFast(FullDBM(mev.ds.update(f)(mev.dec.pure(dimToVarCount(d), ifield.infinity)), mev))

      def flipVar[S <: DBMState, A](vi: VarIndex)(m: CFastDBM[M, SM, S, A])
                                   (implicit ifield: InfField[A]): CFastDBM[M, SM, S, A] =
        Utils.mapFastDBM[M, SM, S, A](Utils.mapInnerMatrix(mev.ds.flipVar(vi)))(m)

      def dbmUnion[S <: DBMState, A](m1: CFastDBM[M, SM, S, A], m2: CFastDBM[M, SM, S, A])
                                    (implicit ifield: InfField[A]): CFastDBM[M, SM, S, A] = {

        def aux(m1: FastDBM[M, SM, A], m2: FastDBM[M, SM, A]): FastDBM[M, SM, A] =
          (m1, m2) match {
            case (FullDBM(dbm1, mev), FullDBM(dbm2, _)) =>
              FullDBM(mev.ds.dbmUnion(dbm1, dbm2), mev)
            case (DecomposedDBM(dbm1, comps1, mev), DecomposedDBM(dbm2, comps2, _)) =>
              val vars1 = comps1.foldRight(Seq[VarIndex]())(_ ++ _).toSet
              val vars2 = comps2.foldRight(Seq[VarIndex]())(_ ++ _).toSet
              val vars = vars1 intersect vars2
              val newComps = comps1.map(_.filter(vars.contains(_)))
              val matrices = newComps.map(c => {
                  val sub1 = mev.dec.extract(c)(dbm1)
                  val sub2 = mev.dec.extract(c)(dbm2)
                  mev.sub.dbmUnion(sub1, sub2)
                })
              val newMat = matrices.foldLeft(dbm1)((mat, subMat) =>
                  mev.dec.pour(subMat)(mat)
                )
              DecomposedDBM(newMat, newComps, mev)
            case (dbm1 @ DecomposedDBM(_, _, _), dbm2 @ FullDBM(_, _)) =>
              aux(dbm1.toFull, dbm2)
            case (dbm1 @ FullDBM(_, _), dbm2 @ DecomposedDBM(_, _, _)) =>
              aux(dbm1, dbm2.toFull)
          }

        (m1, m2) match {
          case (BottomFast(nOfVars), _) =>
            BottomFast(nOfVars)
          case (_, BottomFast(nOfVars)) =>
            BottomFast(nOfVars)
          case (CFast(dbm1), CFast(dbm2)) =>
            CFast(aux(dbm1, dbm2))
          case (CFast(dbm1), NCFast(dbm2)) =>
            // WTF doesn't work TODO
            // NCFast(aux(dbm1, dbm2))
            dbmUnion(m2, m1)
          case (NCFast(dbm1), CFast(dbm2)) =>
            NCFast(aux(dbm1, dbm2))
          case (NCFast(dbm1), NCFast(dbm2)) =>
            NCFast(aux(dbm1, dbm2))
        }
      }

      def addScalarOnVar[S <: DBMState, A]
        (vi: VarIndex, const: A)(m: CFastDBM[M, SM, S, A])
        (implicit ifield: InfField[A]): CFastDBM[M, SM, S, A] =
        Utils.mapFastDBM[M, SM, S, A](
          Utils.mapInnerMatrix(mev.ds.addScalarOnVar(vi, const)))(m)

      def isBottomDBM[A, S <: DBMState](m: CFastDBM[M, SM, S, A])
                                      (implicit ifield: InfField[A]): Boolean =
        Utils.cfastInnerMatrix(m).isEmpty

      def widening[A, R <: DBMState, S <: DBMState]
        (dbm1: CFastDBM[M, SM, R, A], dbm2: CFastDBM[M, SM, S, A])
        (implicit ifield: InfField[A]): ExistsM[A] = {

        def aux(m1: FastDBM[M, SM, A], m2: FastDBM[M, SM, A]): FastDBM[M, SM, A] =
          (m1, m2) match {
            case (FullDBM(dbm1, mev), FullDBM(dbm2, _)) =>
              FullDBM(mev.ds.widening(dbm1, dbm2), mev)
            case (DecomposedDBM(dbm1, comps1, mev), DecomposedDBM(dbm2, comps2, _)) =>
              val vars1 = comps1.foldRight(Seq[VarIndex]())(_ ++ _).toSet
              val vars2 = comps2.foldRight(Seq[VarIndex]())(_ ++ _).toSet
              val vars = vars1 intersect vars2
              val newComps = comps1.map(_.filter(vars.contains(_)))
              val matrices = newComps.map(c => {
                  val sub1 = mev.dec.extract(c)(dbm1)
                  val sub2 = mev.dec.extract(c)(dbm2)
                  mev.sub.widening(sub1, sub2)
                })
              val newMat = matrices.foldLeft(dbm1)((mat, subMat) =>
                  mev.dec.pour(subMat)(mat)
                )
              DecomposedDBM(newMat, newComps, mev)
            case (dbm1 @ DecomposedDBM(_, _, _), dbm2 @ FullDBM(_, _)) =>
              aux(dbm1.toFull, dbm2)
            case (dbm1 @ FullDBM(_, _), dbm2 @ DecomposedDBM(_, _, _)) =>
              aux(dbm1, dbm2.toFull)
          }

        (dbm1, dbm2) match {
          case (BottomFast(nOfVars), _) =>
            Utils.packEx(BottomFast(nOfVars))
          case (_, BottomFast(nOfVars)) =>
            Utils.packEx(BottomFast(nOfVars))
          case (CFast(m1), CFast(m2)) =>
            Utils.packEx(NCFast(aux(m1, m2)))
          case (CFast(m1), NCFast(m2)) =>
            Utils.packEx(NCFast(aux(m1, m2)))
          case (NCFast(m1), CFast(m2)) =>
            Utils.packEx(NCFast(aux(m1, m2)))
          case (NCFast(m1), NCFast(m2)) =>
            Utils.packEx(NCFast(aux(m1, m2)))
        }
      }

      def narrowing[A, R <: DBMState, S <: DBMState]
        (dbm1: CFastDBM[M, SM, R, A], dbm2: CFastDBM[M, SM, S, A])
        (implicit ifield: InfField[A]): ExistsM[A] = {
        val nOfVars = Utils.nOfVars(dbm1)
        (Utils.cfastInnerMatrix(dbm1), Utils.cfastInnerMatrix(dbm2)) match {
          case (None, _) => Utils.packEx(BottomFast(nOfVars))
          case (_, None) => Utils.packEx(BottomFast(nOfVars))
          case (Some(m1), Some(m2)) =>
            val newMat = mev.ds.narrowing(m1, m2)
            Utils.packEx(NCFast(FullDBM(newMat, mev)))
        }
      }

      def isTopDBM[A, S <: DBMState](dbm: CFastDBM[M,SM, S,A])
        (implicit ifield: InfField[A]): Boolean = ???

      def addVariable[S <: DBMState, A](dbm: CFastDBM[M,SM, S,A])
          (implicit ifield: InfField[A]): CFastDBM[M,SM, S,A] = dbm match {
            case BottomFast(n) => BottomFast(addOne(n))
            case _ =>
              Utils.mapFastDBM[M, SM, S, A](
                Utils.mapInnerMatrix[M, SM, A](mev.dec.addVariable))(dbm)
          }

      def decideState[S <: DBMState, A](dbm: CFastDBM[M,SM, S,A]):
          DBMIxed[({ type T[W,B] = CFastDBM[M, SM, W, B]})#T, A] =  dbm match {
          case m @ BottomFast(_) =>
            CIxed[({ type T[W,B] = CFastDBM[M, SM, W, B]})#T, A](m)
          case m @ CFast(_) =>
            CIxed[({ type T[W,B] = CFastDBM[M, SM, W, B]})#T, A](m)
          case m @ NCFast(_) =>
            NCIxed[({ type T[W,B] = CFastDBM[M, SM, W, B]})#T, A](m)
        }

      def deleteVariable[S <: DBMState, A](v: VarIndex)(dbm: CFastDBM[M,SM,S,A])
        (implicit ifield: InfField[A]): CFastDBM[M,SM,S,A] = {
        def f(w: VarIndex): Option[VarIndex] =
          if (w == v) None else
            if (w < v) Some(w) else
              Some(VarIndex(w.i - 1))
        mapVariables(f)(dbm)
      }

      def mapVariables[S <: DBMState, A]
        (f: VarIndex => Option[VarIndex])(dbm: CFastDBM[M,SM,S,A])
        (implicit ifield: InfField[A]): CFastDBM[M,SM,S,A] = dbm match {
        case BottomFast(n) =>
          val newN = allVars(n).count(f(_).isDefined)
          BottomFast(VarCount(newN))
        case _ =>
          Utils.mapFastDBM[M, SM, S, A](fast =>
            Utils.mapInnerMatrix[M, SM, A](mev.dec.mapVariables(f))(fast.toFull))(dbm)
      }

      def compare[A](x: ExistsM[A], y: ExistsM[A])
                    (implicit evidence: InfField[A]): Option[Ordering] = {
        val dbm1: Option[M[A]] = Utils.cfastInnerMatrix(x.elem)
        val dbm2: Option[M[A]] = Utils.cfastInnerMatrix(y.elem)
        (dbm1, dbm2) match {
          case (None, None) => Some(EQ)
          case (Some(_), None) => Some(GT)
          case (None, Some(_)) => Some(LT)
          case (Some(m1), Some(m2)) => mev.ds.compare(m1, m2)
        }
      }
    }
}

// ADT of "closable" DBMs in their fast implementation from Vechev et al.
// They are "closable" in the sense that they augment the ADT of fast DBMs with
// the type-level capability of being indexed by their strong closure state.
sealed trait CFastDBM[M[_], SM[_], _, A]
// Constructor of *closed* fast DBMs.
case class CFast[M[_], SM[_], A](m: FastDBM[M, SM, A]) extends CFastDBM[M, SM, Closed, A]
// Constructor of *non-closed* fast DBMs.
case class NCFast[M[_], SM[_], A](m: FastDBM[M, SM, A]) extends CFastDBM[M, SM, NonClosed, A]
case class BottomFast[M[_], SM[_], A](nOfVars: VarCount) extends CFastDBM[M, SM, Closed, A]

object Utils {

  // TODO: find a better definition.
  // For example, consider nuff decomposed those components
  // with size <= 50% of the total.
  def nuffDecomposed(indepComps: Seq[Seq[VarIndex]], count: VarCount): Boolean =
    indepComps.length >= 2

  def packEx[M[_], SM[_], S <: DBMState, A](fastDBM: CFastDBM[M, SM, S, A])
  : ExistsDBM[({ type T[W] = CFastDBM[M, SM, W, A]})#T] =
    MkEx[S, ({ type T[S] = CFastDBM[M, SM, S, A]})#T](fastDBM)

  def nOfVars[M[_], SM[_], S, A](dbm: CFastDBM[M, SM, S, A])
                         (implicit mev: MEvidence[M, SM]): VarCount =
    dbm match {
      case CFast(m: FastDBM[M, SM, A]) => mev.ds.nOfVars(fastInnerMatrix(m))
      case NCFast(m: FastDBM[M, SM, A]) => mev.ds.nOfVars(fastInnerMatrix(m))
      case BottomFast(n) => n
    }

  def fastInnerMatrix[M[_], SM[_], S <: DBMState, A](fdbm: FastDBM[M, SM, A]): M[A] =
    fdbm match {
      case FullDBM(m: M[A], _) => m
      case DecomposedDBM(m: M[A], _, _) => m
    }

  def cfastInnerMatrix[M[_], SM[_], S <: DBMState, A]
    (cfdbm: CFastDBM[M, SM, S, A])
    (implicit mev: MEvidence[M, SM], ifield: InfField[A]): Option[M[A]] = {

    cfdbm match {
      case CFast(m: FastDBM[M, SM, A]) => Some(fastInnerMatrix(m))
      case NCFast(m: FastDBM[M, SM, A]) => Some(fastInnerMatrix(m))
      case BottomFast(_) => None
    }
  }

  // Applying an arbitrary map destroys every knowledge about relations between
  // variables, so we just return DBM with unknown closure state (because it
  // could be either non-closed or bottom, which is closed).
  // WARNING: input and output matrices of f as assumed to be of the same
  // dimension!
  def liftFromInner[M[_], SM[_], S <: DBMState, A]
    (f : M[A] => M[A])(dbm: CFastDBM[M, SM, S, A])
    (implicit mev: MEvidence[M, SM], ifield: InfField[A])
  : ExistsDBM[({ type T[W] = CFastDBM[M, SM, W, A]})#T] =
    Utils.cfastInnerMatrix(dbm).map((inner) =>
      NCFast(FullDBM(f(inner), mev))) match {
      case Some(cm) => packEx(cm)
      case None => packEx(BottomFast(nOfVars(dbm)))
    }

  def mapFastDBM[M[_], SM[_], S <: DBMState, A](f: FastDBM[M, SM, A] => FastDBM[M, SM, A])
    (cfdbm: CFastDBM[M, SM, S, A])(implicit mev: MEvidence[M, SM], ifield: InfField[A]):
      CFastDBM[M, SM, S, A] = {

    cfdbm match {
      case CFast(m: FastDBM[M, SM, A]) => CFast(f(m))
      case NCFast(m: FastDBM[M, SM, A]) => NCFast(f(m))
      case BottomFast(n) => BottomFast(n)
    }
  }

  def mapInnerMatrix[M[_], SM[_], A](f: M[A] => M[A])(dbm: FastDBM[M, SM, A])
                             (implicit mev: MEvidence[M, SM], ifield: InfField[A]): FastDBM[M, SM, A] = {
    dbm match {
      case FullDBM(m, _) => FullDBM(f(m), mev)
      case DecomposedDBM(m, comps, _) => DecomposedDBM(f(m), comps, mev)
    }
  }

}

object FastDbmUtils {

  def calculateComponents[M[_], A](m: M[A], ds: DenseSparse[M])
    (implicit ifield: InfField[A]): List[List[VarIndex]] = {

    import VarIndexOps._, CountOps._

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

sealed trait FastDBM[M[_], SM[_], A] {
  def strongClosure(implicit mev: MEvidence[M, SM], ifield: InfField[A]):
      CFastDBM[M, SM, Closed, A]

  def incrementalClosure(v: VarIndex)(implicit ifield: InfField[A]):
      CFastDBM[M, SM, Closed, A]

  def toFull: FullDBM[M, SM, A]
}

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
      FastDbmUtils.calculateComponents(dbm, mev.ds)

    if (Utils.nuffDecomposed(indepComponents, mev.ds.nOfVars(dbm)))
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
      FastDbmUtils.calculateComponents(completeDBM, mev.ds)

    if (Utils.nuffDecomposed(indepComponents, mev.ds.nOfVars(completeDBM)))
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

}
