package it.unich.jandom.domains.numerical.octagon
// import scalaz.{Applicative, Monoid}

import scala.language.higherKinds

sealed trait DBMState
case class Closed() extends DBMState
case class NonClosed() extends DBMState

// Distinguish integers used as variable indices
case class VarIndex(i: Int)

// Trait of Difference Bound Matrices, indexed by the closure state
// (closed/non-closed) and the type of the elements.
// Most operators require the type of elements to be a ring.
trait DifferenceBoundMatrix[M[_, _]] {

  def get[A](i: Int, j: Int)(m: M[DBMState, A]): Option[A]
  def update[A](f: (Int, Int) => A)(m: M[DBMState, A]): M[DBMState, A]

  // strong closure and incremental closure are assumed to test for emptiness,
  // and return the bottom element in the positive case.
  def strongClosure[A](m: M[DBMState, A])
    (implicit evidence: InfField[A]): M[Closed, A]
  def incrementalClosure[A](v: VarIndex)(m: M[DBMState, A])
    (implicit evidence: InfField[A]): M[Closed, A]
  def bottomDBM[A]: M[Closed, A]
  def isBottomDBM[A](m: M[DBMState, A]): Boolean
  def topDBM[A]: M[Closed, A]

  // dbm union preserves strong closure
  def dbmUnion[S <: DBMState, A](m1: M[S, A], m2: M[S, A])
    (implicit ifield: InfField[A]): M[S, A]

  // dbm intersection is exact regardless of the closure state of the inputs,
  // and it seldomly produces a strongly closed result.
  def dbmIntersection[A](m1: M[DBMState, A], m2: M[DBMState, A])
      (implicit ifield: InfField[A]): M[DBMState, A]

  /////////////////////////////////////////////////////////////////////////////

  // I feel it should be possible to come up with a more elementary and
  // fine-grained closure state-preserving operation in terms of which all of
  // these closure state-preserving operations can be encoded externally.

  // flip a variable, i.e. interpret v := - v
  // this operation preserves the closure state
  def flipVar[S <: DBMState, A](v: VarIndex)
    (m: M[S, A])(implicit ifield: InfField[A]): M[S, A]

  // add scalar on a variable, i.e. interpret v := v + c
  // this operation preserves the closure state
  def addScalarOnVar[S <: DBMState, A](v: VarIndex, c: A)
    (m: M[S, A])(implicit ifield: InfField[A]): M[S, A]

  // forget operator preserves strong closure
  def forget[S <: DBMState, A](v: VarIndex)(m: M[S, A]): M[S, A]

  //////////////////////////////////////////////////////////////////////////////

  def nOfVars[A](m: M[DBMState, A]): Int
}
