package it.unich.jandom.domains.numerical.octagon.fast

import it.unich.jandom.domains.numerical.octagon._

import scala.language.higherKinds
import variables.{VarIndex, VarCount}

case class NNI(nni: Int)

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

trait DenseSparse[M[_]] {

  def get[A](i: Int, j: Int)(m: M[A]): A

  def update[A](f: (Int, Int) => A)(m: M[A]): M[A]
  def update[A](i: Int, j: Int, x: A)(m: M[A]): M[A]

  def dbmUnion[A](m1: M[A], m2: M[A])(implicit e: InfField[A]): M[A]
  def dbmIntersection[A](m1: M[A], m2: M[A])(implicit e: InfField[A]): M[A]

  def widening[A](m1: M[A], m2: M[A])(implicit e: InfField[A]): M[A]
  def narrowing[A](m1: M[A], m2: M[A])(implicit e: InfField[A]): M[A]

  // closure (incremental or not) may be dense or sparse
  def strongClosure[A](m: M[A])(implicit e: InfField[A]): Option[M[A]]
  def incrementalClosure[A](v: VarIndex)(m: M[A])(implicit e: InfField[A]): Option[M[A]]

  def forget[A](v: VarIndex)(m: M[A])(implicit e: InfField[A]): M[A]

  def flipVar[A](v: VarIndex)(m: M[A]): M[A]

  def addScalarOnVar[A](v: VarIndex, c: A)(m: M[A])
                       (implicit ifield: InfField[A]): M[A]

  def compare[A](m1: M[A], m2: M[A])(implicit ifield: InfField[A]): Option[Ordering]

  //////////////////////////////////////////////////////////////////////////////

  def nOfVars[A](m: M[A]): VarCount

  // // Utility to combine two matrices only on a subset of elements
  // def combine[A, B, C](f: (A, B) => C,
  //   xScope: Seq[Int], yScope: Seq[Int])(ma: M[A], mb: M[B]): M[C]

}

trait SubMatrix[M[_]] extends DenseSparse[M] {
  // The DBM must be aware of the variables that has under control.
  def varIndices[A](m: M[A]): Seq[VarIndex]

  override def nOfVars[A](m: M[A]): VarCount = VarCount(varIndices(m).length)
}

trait Decomposable[M[_], DS[_]] extends SubMatrix[DS] {

  def addVariable[A](m: M[A])(implicit ifield: InfField[A]): M[A]
  def deleteVariable[A](m: M[A])(implicit ifield: InfField[A]): M[A]
  def mapVariables[A](f: VarIndex => Option[VarIndex])(m: M[A])
                     (implicit ifield: InfField[A]): M[A]

  def extract[A](is: Seq[VarIndex])(m: M[A]): DS[A]
  def pour[A](source: DS[A])(dest: M[A]): M[A]

  def pure[A](d: VarCount, x: A): M[A]

}
