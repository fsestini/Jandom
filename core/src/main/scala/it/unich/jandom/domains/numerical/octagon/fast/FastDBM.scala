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
  def instance[M[_]] =
    new DifferenceBoundMatrix[({ type T[S, A] = CFastDBM[M, S, A] })#T] {

      def update[S <: DBMState, A](f: (Int, Int) => A)(m: CFastDBM[M, S, A]): ExistsM[A] =
        ???

      def incrementalClosure[S <: DBMState, A](v: VarIndex)
                               (dbm: CFastDBM[M, S, A])
        (implicit evidence: InfField[A]): CFastDBM[M, Closed, A] =
        ???

      def strongClosure[S <: DBMState, A](dbm: CFastDBM[M, S, A])
                          (implicit evidence: InfField[A]): CFastDBM[M, Closed, A] =
        ???

      def forget[S <: DBMState, A](vi: VarIndex)(m: CFastDBM[M, S, A])
                                  (implicit ifield: InfField[A]): CFastDBM[M, S, A] = 
        ???

      def nOfVars[S <: DBMState, A](m: CFastDBM[M, S, A]): Int = ???

      def get[S <: DBMState, A](i: Int, j: Int)(m: CFastDBM[M, S, A]): Option[A] =
        ???

      def dbmIntersection[A, R <: DBMState, S <: DBMState]
        (m1: CFastDBM[M, R, A], m2: CFastDBM[M, S, A])
        (implicit ifield: InfField[A]): ExistsM[A] = 
        ???

      def topDBM[A](nOfVars: Int)(implicit ifield: InfField[A]): CFastDBM[M, Closed, A] = ???
      def bottomDBM[A](nOfVars: Int)(implicit ifield: InfField[A]): CFastDBM[M, Closed, A] = ???
      def fromFun[A](d: Int, f: ((Int, Int) => A))(implicit ifield: InfField[A]): CFastDBM[M, Closed, A] =
        ???
      def flipVar[S <: DBMState, A](vi: VarIndex)(dbm: CFastDBM[M, S, A])
                                   (implicit ifield: InfField[A]): CFastDBM[M, S, A] =
        ???

      def dbmUnion[S <: DBMState, A](m1: CFastDBM[M, S, A], m2: CFastDBM[M, S, A])
                                    (implicit ifield: InfField[A]): CFastDBM[M, S, A] =
        ???

      def addScalarOnVar[S <: DBMState, A](vi: VarIndex, const: A)
                                          (fundbm: CFastDBM[M, S, A])
                                          (implicit ifield: InfField[A]): CFastDBM[M, S, A] = 
        ???

      def isBottomDBM[A, S <: DBMState](m: CFastDBM[M, S, A]): Boolean =
        ???

      def widening[A, R <: DBMState, S <: DBMState]
        (dbm1: CFastDBM[M, R, A], dbm2: CFastDBM[M, S, A])
        (implicit ifield: InfField[A]): ExistsM[A] =
        ???

      def narrowing[A, R <: DBMState, S <: DBMState]
        (dbm1: CFastDBM[M, R, A], dbm2: CFastDBM[M, S, A])
        (implicit ifield: InfField[A]): ExistsM[A] =
        ???
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

case class BottomDBM[M[_], A]() extends FastDBM[M, A]
