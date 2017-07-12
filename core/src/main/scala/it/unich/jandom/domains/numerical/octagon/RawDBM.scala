package it.unich.jandom.domains.numerical.octagon

import scala.language.higherKinds

// Pack all these operations in a trait, to allow for multiple implementations
// (for ex., full vs. Apron-style half matrices, or parallel stuff, or w/e)
trait RawDBM[M[_]] extends Matrix[M] {

  // All closures return the number of non-infinite values in the matrix,
  // and the independent components of the matrix.
  type VarIndex = Int // Maybe change?
  type ClosureRes[A] = (M[A], Int, List[List[VarIndex]])

  // The list on indices is passed when using decomposed DBMs.

  def nOfVars: Int

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
