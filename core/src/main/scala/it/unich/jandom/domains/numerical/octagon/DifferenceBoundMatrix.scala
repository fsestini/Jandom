package it.unich.jandom.domains.numerical.octagon
// import scalaz.{Applicative, Monoid}

import breeze.math.Field

import scala.language.higherKinds
import it.unich.jandom.domains.numerical._

sealed trait DBMState
sealed trait Closed extends DBMState
sealed trait NonClosed extends DBMState

sealed trait ExistsDBM[M[_]] {
  type State <: DBMState
  val elem: M[State]
}

final case class MkEx[S <: DBMState, M[_]](elem: M[S])
  extends ExistsDBM[M] { type State = S }

// Distinguish integers used as variable indices
case class VarIndex(i: Int)

object VarIndexOps {
  sealed trait OctaVarCoeff
  object Positive extends OctaVarCoeff { }
  object Negative extends OctaVarCoeff { }

  def varPlus(v: VarIndex): Int = 2 * v.i
  def varMinus(v: VarIndex): Int = 2 * v.i + 1
  def signed(i: Int): Int = if (i % 2 == 0) i + 1 else i - 1
}

object VarIndexUtils {
  def forSomeVar(
    vars: Seq[VarIndex])(p: VarIndex => Boolean): Option[VarIndex] =
    (vars.map(x => if (p(x)) Some(x) else None).toList).flatten.headOption
  // Evaluation of linear assignment using interval arithmetics.
  def lfAsInterval(v: VarIndex, lf: LinearForm): (Double, Double) = ???
}

// Trait of Difference Bound Matrices, indexed by the closure state
// (closed/non-closed) and the type of the elements.
// Most operators require the type of elements to be a ring.
trait DifferenceBoundMatrix[M[_, _]] {

  type ExistsM[A] = ExistsDBM[({ type T[S] = M[S, A]})#T]

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
  def bottomDBM[A](nOfVars : Int)(implicit ifield: InfField[A]) : M[Closed, A]
  def isBottomDBM[A, S <: DBMState](dbm: M[S, A]): Boolean
  def topDBM[A](nOfVars : Int)(implicit ifield: InfField[A]) : M[Closed, A]
  def fromFun[A] (d: Int, f: ((Int, Int) => A))(implicit ifield: InfField[A]) : M[Closed, A]

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

  def nOfVars[S <: DBMState, A](m: M[S, A]): Int
}
