package it.unich.jandom.domains.numerical.octagon.fast

import breeze.numerics.pow
import it.unich.jandom.domains.numerical.octagon._

import scala.language.higherKinds
import scalaz.{Applicative, Apply, Monoid, Traverse}
import scalaz.std.option._
import scalaz.std.list._

// Sparsity D = 1 - (nni / (2n^2 + 2n))
// Switching between DBMs: the sparsity can increase, for instance during
// widening. Recovering sparsity information and  independent components has a
// quadratic worst case complexity, so we only perform it by piggybacking on
// the closure operator. We also use closure computations as switching points.

object CFDBMInstance {
  def instance[M[_]](implicit ds: DenseSparseDBM[M]) =
    new DifferenceBoundMatrix[({ type T[S, A] = CFastDBM[M, S, A] })#T] {

      def update[S <: DBMState, A](f: (Int, Int) => A)(m: CFastDBM[M, S, A])
                                  (implicit ifield: InfField[A]): ExistsM[A] =
        Utils.liftFromInner(ds.update(f))(m)

      def incrementalClosure[S <: DBMState, A](v: VarIndex)
                               (dbm: CFastDBM[M, S, A])
        (implicit evidence: InfField[A]): CFastDBM[M, Closed, A] =
        ???

      def strongClosure[S <: DBMState, A](dbm: CFastDBM[M, S, A])
                          (implicit evidence: InfField[A]): CFastDBM[M, Closed, A] =
        ???

      def forget[S <: DBMState, A](vi: VarIndex)(m: CFastDBM[M, S, A])
                                  (implicit ifield: InfField[A]): CFastDBM[M, S, A] =
        Utils.mapFastDBM[M, S, A](fast =>
          Utils.mapInnerMatrix[M, A](inner =>
            ds.forget(vi)(inner)
          )(fast)
        )(m)

      def nOfVars[S <: DBMState, A](m: CFastDBM[M, S, A]): Int = nOfVars(m)

      def get[S <: DBMState, A](i: Int, j: Int)(m: CFastDBM[M, S, A])
                               (implicit ifield: InfField[A]): Option[A] =
        Utils.cfastInnerMatrix(m).flatMap((inner) => ds.get(i, j)(inner))

      def dbmIntersection[A, R <: DBMState, S <: DBMState]
        (m1: CFastDBM[M, R, A], m2: CFastDBM[M, S, A])
        (implicit ifield: InfField[A]): ExistsM[A] =
        ???

      def topDBM[A](nOfVars: Int)(implicit ifield: InfField[A]): CFastDBM[M, Closed, A] =
        TopFast(nOfVars)

      def bottomDBM[A](nOfVars: Int)(implicit ifield: InfField[A]): CFastDBM[M, Closed, A] =
        BottomFast(nOfVars)

      def fromFun[A](d: Int, f: ((Int, Int) => A))(implicit ifield: InfField[A]): CFastDBM[M, Closed, A] =
        ???
      def flipVar[S <: DBMState, A](vi: VarIndex)(m: CFastDBM[M, S, A])
                                   (implicit ifield: InfField[A]): CFastDBM[M, S, A] = {

        def aux(m: FastDBM[M, A]): FastDBM[M, A] = m match {
          case FullDBM(dbm, dsdbm) =>
            FullDBM(dsdbm.flipVar(vi)(dbm), dsdbm)
          case DecomposedDBM(completeDBM, indepComponents, rdbm) =>
            val comp = indepComponents.find(c => {
                c.contains(vi)
              })
            comp match {
              case Some(c) =>
                val subMat = rdbm.extract(c)(completeDBM)
                val updated = rdbm.flipVar(vi)(subMat)
                val newMat = rdbm.pour(subMat)(completeDBM)
                DecomposedDBM(newMat, indepComponents, rdbm)
              case None =>
                m
            }
        }

        m match {
          case BottomFast(nOfVars) =>
            m
          case TopFast(nOfVars) =>
            m
          case CFast(dbm) =>
            CFast(aux(dbm))
          case NCFast(dbm) =>
            NCFast(aux(dbm))
        }
      }

      def dbmUnion[S <: DBMState, A](m1: CFastDBM[M, S, A], m2: CFastDBM[M, S, A])
                                    (implicit ifield: InfField[A]): CFastDBM[M, S, A] =
        ???

      def addScalarOnVar[S <: DBMState, A](vi: VarIndex, const: A)
                                          (m: CFastDBM[M, S, A])
                                          (implicit ifield: InfField[A]): CFastDBM[M, S, A] = {

        def aux(m: FastDBM[M, A]): FastDBM[M, A] = m match {
          case FullDBM(dbm, dsdbm) =>
            FullDBM(dsdbm.addScalarOnVar(vi, const)(dbm), dsdbm)
          case DecomposedDBM(completeDBM, indepComponents, rdbm) =>
            val comp = indepComponents.find(c => {
                c.contains(vi)
              })
            comp match {
              case Some(c) =>
                val subMat = rdbm.extract(c)(completeDBM)
                val updated = rdbm.addScalarOnVar(vi, const)(subMat)
                val newMat = rdbm.pour(subMat)(completeDBM)
                DecomposedDBM(newMat, indepComponents, rdbm)
              case None =>
                m
            }
        }

        m match {
          case BottomFast(nOfVars) =>
            m
          case TopFast(nOfVars) =>
            m
          case CFast(dbm) =>
            CFast(aux(dbm))
          case NCFast(dbm) =>
            NCFast(aux(dbm))
        }
      }

      def isBottomDBM[A, S <: DBMState](m: CFastDBM[M, S, A])
                                      (implicit ifield: InfField[A]): Boolean =
        Utils.cfastInnerMatrix(m).isEmpty

      def widening[A, R <: DBMState, S <: DBMState]
        (dbm1: CFastDBM[M, R, A], dbm2: CFastDBM[M, S, A])
        (implicit ifield: InfField[A]): ExistsM[A] =
        ???

      def narrowing[A, R <: DBMState, S <: DBMState]
        (dbm1: CFastDBM[M, R, A], dbm2: CFastDBM[M, S, A])
        (implicit ifield: InfField[A]): ExistsM[A] =
        ???

      def isTopDBM[A, S <: DBMState](dbm: CFastDBM[M,S,A])(implicit ifield: InfField[A]): Boolean =
        dbm match {
          case TopFast(_) => true
          case _ => false
        }

      def addVariable[S <: DBMState, A](dbm: CFastDBM[M,S,A])(implicit ifield: InfField[A]): CFastDBM[M,S,A] = ???
      def decideState[S <: DBMState, A](dbm: CFastDBM[M,S,A]): DBMIxed[({ type T[W,B] = CFastDBM[M, W, B]})#T, A] = ???
      def deleteVariable[S <: DBMState, A](v: VarIndex)(dbm: CFastDBM[M,S,A])(implicit ifield: InfField[A]): CFastDBM[M,S,A] = ???
      def mapVariables[S <: DBMState, A](f: VarIndex => Option[VarIndex])(dbm: CFastDBM[M,S,A])(implicit ifield: InfField[A]): CFastDBM[M,S,A] = ???
      def compare[A](x: ExistsM[A], y: ExistsM[A])(implicit pc: this.PosetConstraint[A]): Option[Ordering] = ???
    }
}

// ADT of "closable" DBMs in their fast implementation from Vechev et al.
// They are "closable" in the sense that they augment the ADT of fast DBMs with
// the type-level capability of being indexed by their strong closure state.
sealed trait CFastDBM[M[_], _, A]
// Constructor of *closed* fast DBMs.
case class CFast[M[_], A](m: FastDBM[M, A]) extends CFastDBM[M, Closed, A]
// Constructor of *non-closed* fast DBMs.
case class NCFast[M[_], A](m: FastDBM[M, A]) extends CFastDBM[M, NonClosed, A]
case class TopFast[M[_], A](nOfVars: Int) extends CFastDBM[M, Closed, A]
case class BottomFast[M[_], A](nOfVars: Int) extends CFastDBM[M, Closed, A]

object Utils {

  def packEx[M[_], S <: DBMState, A](fastDBM: CFastDBM[M, S, A])
  : ExistsDBM[({ type T[W] = CFastDBM[M, W, A]})#T] =
    MkEx[S, ({ type T[S] = CFastDBM[M, S, A]})#T](fastDBM)

  def nOfVars[M[_], S, A](dbm: CFastDBM[M, S, A])
                         (implicit ds: DenseSparseDBM[M]): Int =
    dbm match {
      case CFast(m: FastDBM[M, A]) => ds.nOfVars(fastInnerMatrix(m))
      case NCFast(m: FastDBM[M, A]) => ds.nOfVars(fastInnerMatrix(m))
      case TopFast(n) => n
      case BottomFast(n) => n
    }

  def fastInnerMatrix[M[_], S <: DBMState, A](fdbm: FastDBM[M, A]): M[A] =
    fdbm match {
      case FullDBM(m: M[A], _) => m
      case DecomposedDBM(m: M[A], _, _) => m
    }

  def cfastInnerMatrix[M[_], S <: DBMState, A]
    (cfdbm: CFastDBM[M, S, A])(implicit ds: DenseSparseDBM[M], ifield: InfField[A]): Option[M[A]] = {

    cfdbm match {
      case CFast(m: FastDBM[M, A]) => Some(fastInnerMatrix(m))
      case NCFast(m: FastDBM[M, A]) => Some(fastInnerMatrix(m))
      case TopFast(n) => Some(ds.pure(ifield.infinity))
      case BottomFast(_) => None
    }
  }

  // Applying an arbitrary map destroys every knowledge about relations between
  // variables, so we just return DBM with unknown closure state (because it
  // could be either non-closed or bottom, which is closed).
  // WARNING: input and output matrices of f as assumed to be of the same
  // dimension!
  def liftFromInner[M[_], S <: DBMState, A](f : M[A] => M[A])
                                           (dbm: CFastDBM[M, S, A])
                                           (implicit ds: DenseSparseDBM[M],
                                            ifield: InfField[A])
  : ExistsDBM[({ type T[W] = CFastDBM[M, W, A]})#T] =
    Utils.cfastInnerMatrix(dbm).map((inner) =>
      NCFast(FullDBM(f(inner), ds))) match {
      case Some(cm) => packEx(cm)
      case None => packEx(BottomFast(nOfVars(dbm)))
    }

  def mapFastDBM[M[_], S <: DBMState, A](f: FastDBM[M, A] => FastDBM[M, A])
    (cfdbm: CFastDBM[M, S, A])(implicit ds: DenseSparseDBM[M], ifield: InfField[A]): CFastDBM[M, S, A] = {

    cfdbm match {
      case CFast(m: FastDBM[M, A]) => CFast(f(m))
      case NCFast(m: FastDBM[M, A]) => NCFast(f(m))
      case TopFast(n) => CFast(f(FullDBM(ds.pure(ifield.infinity), ds)))
      case BottomFast(n) => BottomFast(n)
    }
  }

  def mapInnerMatrix[M[_], A](f: M[A] => M[A])(dbm: FastDBM[M, A])
                             (implicit ds: DenseSparseDBM[M], ifield: InfField[A]): FastDBM[M, A] = {
    dbm match {
      case FullDBM(m, _) => FullDBM(f(m), ds)
      case DecomposedDBM(m, comps, _) => DecomposedDBM(f(m), comps, ds)
    }
  }

}

object Lol {

  def nuffDecomposed(is: List[List[VarIndex]]): Boolean = ???
  def nuffSparse(is: NNI): Boolean = ???

}

sealed trait FastDBM[M[_], A]

// Full DBMs are fast DBMs that are not decomposed, i.e., they can be either
// dense or sparse.
// Sparsity details, including when to switch between dense and sparse
// representation, is supposed to be handled by the specific implementation of
// the the DenseSparse trait/typeclass.
// An even better thing to do (time permitting) could be to use a suitable
// abstract trait of DBMs that does not talk about sparsity at all (who cares
// if the implementations use a dense/sparse representation anyway, as long as
// they provide the needed DBM-like operations?)
case class FullDBM[M[_], A](dbm: M[A], dsdbm: DenseSparseDBM[M]) extends FastDBM[M, A]

// We store the independent components as a linked list of linked lists of
// variable indices.
// Independent components correspond to submatrices in the complete DBM,
// and these can be dense or sparse. The exact sparsity is computed
// on-the-fly before closure.

// Octagon operators on the decomposed type:
// Apart from closure (which is specialized according to sparsity), we apply the
// standard operators independently on each submatrix.

// NNI for decomposed DBMs is computed as the sum of its submatrices' NNI.
case class DecomposedDBM[M[_], A](completeDBM: M[A],
                                  indepComponents: Seq[Seq[VarIndex]],
                                  rdbm: DenseSparseDBM[M]) extends FastDBM[M, A] {

  // Skeleton of strong closure for decomposed matrices.
  def decStrongClosure(m: DecomposedDBM[M, A])
                      (implicit ifield: InfField[A])
  : (CFastDBM[M, Closed, A], NNI, List[List[VarIndex]]) = {

    val submatrices = indepComponents.map(seq => rdbm.extract(seq)(completeDBM))
    Applicative[Option].sequence(
      submatrices.map(m => rdbm.strongClosure(m)).toList) match {

      case Some(closedSubs) => {
        val (newMatrix, newNNI, newIndepComps) =
          closedSubs.foldRight((completeDBM, 0, indepComponents))(
            (x, y) => (x, y) match {
              case ((sub, NNI(nni), _), (full, fullNNI, is)) =>
                (rdbm.pour(sub)(full), nni + fullNNI, is)})
        val actualNewIndepComponents = ???
        val m =
          if (Lol.nuffDecomposed(actualNewIndepComponents))
            CFast(DecomposedDBM(newMatrix, actualNewIndepComponents, rdbm))
          else
            // might be that we have to compute the index of non-infinite terms
            // during the strong closure above, and pass them to the dbm
            // constructor.
            CFast(FullDBM(newMatrix, rdbm))

        (m, NNI(newNNI), actualNewIndepComponents)
      }
      case None => (BottomFast(rdbm.varIndices(completeDBM).length), ???, Nil)
    }
  }
}
