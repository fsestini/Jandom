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
  def instance[M[_]]: DifferenceBoundMatrix[
    ({ type T[S, A] = CFastDBM[M, S, A] })#T] = ???
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
case class DenseDBM[M[_], A](m: M[A], rdbm: DenseSparseDBM[M]) extends FastDBM[M, A] { }
case class SparseDBM[M[_], A](m: M[A], rdbm: DenseSparseDBM[M]) extends FastDBM[M, A] { }

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
    (implicit ifield: InfField[A]):
      Option[(CFastDBM[M, Closed, A], NNI, List[List[VarIndex]])] = {
    val subs: Seq[M[A]] =
      indepComponents.map(seq => rdbm.extract(seq)(completeDBM))
    val closedSubs: List[Option[(M[A], NNI, List[List[VarIndex]])]] =
      subs.map(m => rdbm.strongClosure(m)).toList
    val actualClosedSubs: Option[List[(M[A], NNI, List[List[VarIndex]])]] =
      Applicative[Option].sequence(closedSubs)
    actualClosedSubs match {
      case Some(cloSu) => {
        val (newM, newNNI, newIs) = cloSu.foldRight((completeDBM, 0, indepComponents))(
          (x, y) => (x, y) match {
            case ((sub, NNI(nni), _), (full, fullNNI, is)) =>
              (rdbm.pour(sub)(full), nni + fullNNI, is)})
        // Computation of the new list of independent components is non-trivial.
        val actualNewIndepComponents = ???
        val m: CFastDBM[M, Closed, A] =
          if (Lol.nuffDecomposed(actualNewIndepComponents)) {
            CFast(DecomposedDBM(newM, actualNewIndepComponents, rdbm))
          } else if (Lol.nuffSparse(NNI(newNNI))) {
            // might be that we have to compute the index of non-infinite terms
            // during the strong closure above, and pass them to the
            // sparse dbm constructor.
            CFast(SparseDBM(newM, rdbm))
          } else {
            CFast(DenseDBM(newM, rdbm))
          }
        Some((m, NNI(newNNI), actualNewIndepComponents))
      }
      case None => None
    }
  }
}
