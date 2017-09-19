package it.unich.jandom.domains.numerical.octagon
// import scalaz.{Applicative, Monoid}

import breeze.math.Field

import scala.language.higherKinds
import it.unich.jandom.domains.numerical._
import it.unich.jandom.domains.numerical.octagon.variables.VarIndex
import it.unich.jandom.domains.numerical.octagon.variables.VarCount
import it.unich.jandom.domains.numerical.octagon.variables.Dimension

sealed trait DBMState
sealed trait Closed extends DBMState
sealed trait NonClosed extends DBMState

sealed trait ExistsDBM[M[_]] {
  type State <: DBMState
  val elem: M[State]
}

final case class MkEx[S <: DBMState, M[_]](elem: M[S])
  extends ExistsDBM[M] { type State = S }

sealed trait DBMIxed[M[_,_], A]
case class CIxed[M[_,_], A](x: M[Closed, A]) extends DBMIxed[M, A]
case class NCIxed[M[_,_], A](x: M[NonClosed, A]) extends DBMIxed[M, A]

// Trait of Difference Bound Matrices, indexed by the closure state
// (closed/non-closed) and the type of the elements.
// Most operators require the type of elements to be a ring.
trait DifferenceBoundMatrix[M[_, _]]
  extends Poset1[({ type T[A] = ExistsDBM[({ type Q[S] = M[S, A]})#Q]})#T] {

  type ExistsM[A] = ExistsDBM[({ type T[S] = M[S, A]})#T]

  def decideState[S <: DBMState, A](dbm: M[S, A]): DBMIxed[M, A]

  // Returns None is the DBM is bottom. Otherwise, Some(element).

  def get[S <: DBMState, A](i: Int, j: Int)(m: M[S, A])
                           (implicit ifield: InfField[A]): Option[A]
  def update[S <: DBMState, A](f: (Int, Int) => A)(m: M[S, A])
                              (implicit ifield: InfField[A]): ExistsM[A]

  // strong closure and incremental closure are assumed to test for emptiness,
  // and return the bottom element in the positive case.
  def strongClosure[S <: DBMState, A](m: M[S, A])
    (implicit evidence: InfField[A]): M[Closed, A]
  def incrementalClosure[S <: DBMState, A](v: VarIndex)(m: M[S, A])
    (implicit evidence: InfField[A]): M[Closed, A]
  def widening[A, S <: DBMState, T <: DBMState]
    (m1: M[S, A], m2: M[T, A])(implicit ifield: InfField[A]): ExistsM[A]
  def narrowing[A, S <: DBMState, T <: DBMState]
    (m1: M[S, A], m2: M[T, A])(implicit ifield: InfField[A]): ExistsM[A]
  def bottomDBM[A](nOfVars : VarCount)(implicit ifield: InfField[A]) : M[Closed, A]

  // Tests whether dbm is a top DBM, i.e. if it has all elements == +\infty.
  // It is provided as primitive, and not implemented in terms of get,
  // to allow more efficient implementations.
  def isTopDBM[A, S <: DBMState](dbm: M[S, A])(implicit ifield: InfField[A]): Boolean
  def isBottomDBM[A, S <: DBMState](dbm: M[S, A])(implicit ifield: InfField[A]): Boolean

  // Returns a *closed* version of the top DBM, namely one with +\infty on
  // all elements except the diagonal, which is set to 0.
  // This differs from the *real* top DBM, which has +\infty on all elements.
  // In particular, isTop(topDBM(n)) == false.
  def topDBM[A](nOfVars : VarCount)(implicit ifield: InfField[A]) : M[Closed, A]
  def fromFun[A] (d: Dimension, f: ((Int, Int) => A))(implicit ifield: InfField[A]) : M[Closed, A]

  // dbm union preserves strong closure
  def dbmUnion[S <: DBMState, A](m1: M[S, A], m2: M[S, A])
    (implicit ifield: InfField[A]): M[S, A]

  // dbm intersection is exact regardless of the closure state of the inputs,
  // and it seldomly produces a strongly closed result.
  def dbmIntersection[A, S <: DBMState, T <: DBMState]
  (m1: M[S, A], m2: M[T, A])
  (implicit ifield: InfField[A]): ExistsM[A]

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
  def forget[S <: DBMState, A](v: VarIndex)(m: M[S, A])
                              (implicit ifield: InfField[A]): M[S, A]

  //////////////////////////////////////////////////////////////////////////////

  def nOfVars[S <: DBMState, A](m: M[S, A]): VarCount
  def addVariable[S <: DBMState, A](dbm: M[S, A])(implicit ifield: InfField[A]): M[S, A]
  def deleteVariable[S <: DBMState, A](v: VarIndex)(dbm: M[S, A])(implicit ifield: InfField[A]): M[S, A]
  def mapVariables[S <: DBMState, A](f: VarIndex => Option[VarIndex])
                                    (dbm: M[S, A])
                                    (implicit ifield: InfField[A]): M[S, A]
}
