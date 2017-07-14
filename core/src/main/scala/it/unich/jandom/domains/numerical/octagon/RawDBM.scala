package it.unich.jandom.domains.numerical.octagon

import scala.language.higherKinds

// This trait describes DBM matrices that can be used as building blocks in the
// implementation of the "fast" version of abstract octagons.

// This trait and DifferenceBoundMatrix live at different levels of abstraction,
// and thus should be kept separated. As an example, the ADT encoding DBMs in
// the "fast" version of Vechev et al. is an instance of DifferenceBoundMatrix,
// but may use several RawDBM objects internally. In particular, the Decomposed
// constructor of the ADT is a DifferenceBoundMatrix as a whole, but is
// internally implemented as a collection of RawDBMs.
// Notice that, for example, a strong closure operation or "forget" operation on
// a Decomposed DBM is realized as the independent execution of the same
// operation on the RawDBM objects, so DifferenceBoundMatrix and RawDBM
// may need to share some operations.
// In the future, we may try to collect these common set of operations in a
// common trait that is extended by both.

// A differenza di Diff..., this trait only makes sense in the implementation of
// "fast" octagons (it is aware of the dense/sparse representations).
// Diff... instead is directly accessed by ABstractOctagon, and describes a
// generic DBM encoding octagon constraints.

// Pack all these operations in a trait, to allow for multiple implementations
// (for ex., full vs. Apron-style half matrices, or parallel stuff, or w/e)
trait RawDBM[M[_]] extends Matrix[M] {

  // All closures return the number of non-infinite values in the matrix,
  // and the independent components of the matrix.
  type VarIndex = Int // Maybe change?
  type ClosureRes[A] = (M[A], Int, List[List[VarIndex]])

  // The list on indices is passed when using decomposed DBMs.

  def nOfVars[A](m: M[A]): Int

  def denseStrongClosure[A](m: M[A])(implicit e: InfField[A]): ClosureRes[A]
  def denseStrongClosure[A](m: M[A], indices: Seq[Int])(implicit e: InfField[A]): ClosureRes[A]
  def sparseStrongClosure[A](m: M[A])(implicit e: InfField[A]): ClosureRes[A]
  def sparseStrongClosure[A](m: M[A], indices: Seq[Int])(implicit e: InfField[A]): ClosureRes[A]

  def denseIncrementalClosure[A](m: M[A])(implicit e: InfField[A]): ClosureRes[A]
  def sparseIncrementalClosure[A](m: M[A])(implicit e: InfField[A]): ClosureRes[A]

  // Utility to combine two matrices only on a subset of elements
  def combine[A, B, C](f: (A, B) => C,
    xScope: Seq[Int], yScope: Seq[Int])(ma: M[A], mb: M[B]): M[C]
}
