package it.unich.jandom.domains.numerical.octagon

import scalaz.{Applicative, Apply, Monoid}
import scalaz.std.option._, scalaz.std.list._

/**
  * Created by fsestini on 7/13/17.
  */

// FunMatrix-based raw DBM implementation
sealed trait FunDBM[S, A] {
  def liftFromInner(f: FunMatrix[A] => FunMatrix[A]): FunDBM[S, A]
  def union(other: FunDBM[S, A])(implicit infField: InfField[A]): FunDBM[S, A]
  val innerMatrix: Option[FunMatrix[A]]
  val noOfVariables: Int
}

case class ClosedFunDBM[A](nOfVars: Int, m: FunMatrix[A]) extends FunDBM[Closed, A] {
  override def liftFromInner(f: (FunMatrix[A]) => FunMatrix[A]): FunDBM[Closed, A] =
    ClosedFunDBM(nOfVars, f(m))

  override def union(other: FunDBM[Closed, A])
                    (implicit infField: InfField[A]): FunDBM[Closed, A] = other match {
    case ClosedFunDBM(n2, m2) => {
      val me: Matrix[FunMatrix] = FunMatrixMatrixInstance.funMatrixIsMatrix
      require(nOfVars == n2)
      ClosedFunDBM(nOfVars, me.combine(infField.max)(m, m2))
    }
    case TopFunDBM(n2, infty) => { require(nOfVars == n2) ; TopFunDBM(n2, infty) }
    case BottomFunDBM(n2) => { require(nOfVars == n2) ; this }
  }

  val innerMatrix: Option[FunMatrix[A]] = Some(m)
  val noOfVariables: Int = nOfVars
}

case class NonClosedFunDBM[A](nOfVars: Int, m: FunMatrix[A]) extends FunDBM[NonClosed, A] {
  override def liftFromInner(f: (FunMatrix[A]) => FunMatrix[A]): FunDBM[NonClosed, A] =
    NonClosedFunDBM(nOfVars, f(m))
  def union(other: FunDBM[NonClosed, A])
           (implicit infField: InfField[A]): FunDBM[NonClosed, A] = other match {
    case NonClosedFunDBM(n2, m2) => {
      val me: Matrix[FunMatrix] = FunMatrixMatrixInstance.funMatrixIsMatrix
      require(nOfVars == n2)
      NonClosedFunDBM(nOfVars, me.combine(infField.max)(m, m2))
    }
  }

  val innerMatrix: Option[FunMatrix[A]] = Some(m)
  val noOfVariables: Int = nOfVars
}

case class TopFunDBM[A](nOfVars: Int, infty: A) extends FunDBM[Closed, A] {
  override def liftFromInner(f: (FunMatrix[A]) => FunMatrix[A]): FunDBM[Closed, A] =
    TopFunDBM[A](nOfVars, infty)
  def union(other: FunDBM[Closed, A])(implicit infField: InfField[A]): FunDBM[Closed, A] = this

  val innerMatrix: Option[FunMatrix[A]] =
    Some(FunMatrixMatrixInstance.funMatrixIsMatrix.pure(2 * nOfVars, infty))
  val noOfVariables: Int = nOfVars
}

case class BottomFunDBM[A](nOfVars: Int) extends FunDBM[Closed, A] {
  override def liftFromInner(f: (FunMatrix[A]) => FunMatrix[A]): FunDBM[Closed, A] =
    BottomFunDBM[A](nOfVars)
  def union(other: FunDBM[Closed, A])(implicit infField: InfField[A]): FunDBM[Closed, A] = other

  val innerMatrix: Option[FunMatrix[A]] = None
  val noOfVariables: Int = nOfVars
}

object FunDBMInstance {
  val me: Matrix[FunMatrix] = FunMatrixMatrixInstance.funMatrixIsMatrix

  implicit val funDBM: DifferenceBoundMatrix[FunDBM] = new DifferenceBoundMatrix[FunDBM] {
    def update[A](f: (Int, Int) => A)(m: FunDBM[DBMState, A]): FunDBM[DBMState, A] =
      m.liftFromInner(me.update(f))

    def incrementalClosure[A](v: VarIndex)
                             (m: FunDBM[DBMState, A])
                             (implicit evidence: InfField[A]): FunDBM[Closed, A] = ???

    def strongClosure[A](dbm: FunDBM[DBMState, A])
                        (implicit evidence: InfField[A]): FunDBM[Closed, A] =
      dbm.innerMatrix match {
        case Some(m) =>
          ClosedFunDBM(
            dbm.noOfVariables,
            BagnaraStrongClosure.strongClosure(dbm.noOfVariables)(m))
        case None => BottomFunDBM(dbm.noOfVariables)
      }

    def forget[S <: DBMState, A](vi: VarIndex)(m: FunDBM[S, A])
                                (implicit ifield: InfField[A]): FunDBM[S, A] = {
      val v = vi.i
      m.liftFromInner((inner) => {
        val f: (Int, Int) => A = (i, j) => {
          if (i != 2 * v - 1 && i != 2 * v && j != 2 * v - 1 && j != 2 * v)
            me.get(i, j)(inner)
          else if (i == j && (i == 2 * v - 1 || i == 2 * v))
            ifield.zero
          else ifield.infinity
        }
        me.update(f)(inner)
      })
    }

    def nOfVars[A](m: FunDBM[DBMState, A]): Int = m.noOfVariables

    def get[A](i: Int, j: Int)(m: FunDBM[DBMState, A]): Option[A] =
      m.innerMatrix.map((mmm) => me.get(i, j)(mmm))

    def dbmIntersection[A, S <: DBMState, T <: DBMState]
    (m1: FunDBM[S, A], m2: FunDBM[T, A])
    (implicit ifield: InfField[A])
    : FunDBM[W, A] forSome { type W <: DBMState } = {
      require(m1.noOfVariables == m2.noOfVariables)
      val o = for {
        mm1 <- m1.innerMatrix
        mm2 <- m2.innerMatrix
      } yield NonClosedFunDBM(m1.noOfVariables, me.combine(ifield.min)(mm1, mm2))
      o match {
        case Some(matrix) => matrix
        case None => BottomFunDBM(m1.noOfVariables)
      }
    }

    def topDBM[A](nOfVars: Int)(implicit ifield: InfField[A]): FunDBM[Closed, A] = TopFunDBM(nOfVars, ifield.infinity)
    def bottomDBM[A](nOfVars: Int)(implicit ifield: InfField[A]): FunDBM[Closed, A] = BottomFunDBM(nOfVars)

    private def signed(i: Int): Int = ???

    def flipVar[S <: DBMState, A](vi: VarIndex)(dbm: FunDBM[S, A])
                                 (implicit ifield: InfField[A]): FunDBM[S, A] = {
      val v = vi.i
      dbm.liftFromInner((inner) => {
        val f: (Int, Int) => A = (i, j) => {
          if (i == 2 * v - 1 || i == 2 * v) {
            if (j == 2 * v - 1 || j == 2 * v)
              me.get(signed(i), signed(i))(inner)
            else
              me.get(signed(i), j)(inner)
          } else {
            if (j == 2 * v - 1 || j == 2 * v)
              me.get(i, signed(j))(inner)
            else
              me.get(i, j)(inner)
          }
        }
        me.update(f)(inner)
      })
    }

    def dbmUnion[S <: DBMState, A](m1: FunDBM[S, A], m2: FunDBM[S, A])
                                  (implicit ifield: InfField[A]): FunDBM[S, A] = m1.union(m2)


    def addScalarOnVar[S <: DBMState, A](vi: VarIndex, const: A)
                                        (fundbm: FunDBM[S, A])
                                        (implicit ifield: InfField[A]): FunDBM[S, A] = {
      val v = vi.i
      fundbm.liftFromInner((dbm) => {
        val f: (Int, Int) => A = (i, j) => {
          val g1 = (i == 2 * v - 1 && j != 2 * v - 1 && j != 2 * v) ||
                   (j == 2 * v && i != 2 * v - 1 && i != 2 * v)
          val g2 = (i != 2 * v - 1 && i != 2 * v && j == 2 * v - 1) ||
                   (j != 2 * v - 1 && j != 2 * v && i == 2 * v)
          val g3 = i == 2 * v - 1 && j == 2 * v
          val g4 = i == 2 * v && j == 2 * v - 1
          if (g1) ifield.-(me.get(i, j)(dbm), const) else
            if (g2) ifield.+(me.get(i, j)(dbm), const) else
              if (g3) ifield.-(me.get(i, j)(dbm), ifield.double(const)) else
                if (g4) ifield.+(me.get(i, j)(dbm), ifield.double(const)) else
                  me.get(i, j)(dbm)
        }
        me.update(f)(dbm)
      })
    }

    def isBottomDBM[A, S <: DBMState](m: FunDBM[S, A]): Boolean = m match {
      case BottomFunDBM(_) => true
      case _ => false
    }
  }
}

// This is the simplified strong closure algorithm from
// Bagnata et al., Widening Operators for Weakly-Relational Numeric Abstractions.
//
// It is a classical Floyd-Warshall, followed by a single strong coherence step.
//
// for k = 0 to 2*n - 1
//   for i = 0 to 2*n - 1
//     for j = 0 to 2*n - 1
//       m(i,j) = min( m(i,j), m(i,k) + m(k,j) )
//
// for i = 0 to 2*n - 1
//   for j = 0 to 2*n - 1
//     m(i,j) = min( m(i,j), m(i, signed(i)) + m(signed(j), j) / 2)
object BagnaraStrongClosure {
  private val me: Matrix[FunMatrix] = FunMatrixMatrixInstance.funMatrixIsMatrix

  private def signed(i: Int) = if (i % 2 == 0) i + 1 else i - 1

  def strongClosure[A](nOfVars: Int)(dbm: FunMatrix[A])
                      (implicit ifield: InfField[A]): FunMatrix[A] = {
    var x = dbm
    for (k <- 0 to 2 * nOfVars - 1)
      for (i <- 0 to 2 * nOfVars - 1)
        for (j <- 0 to 2 * nOfVars - 1) {
          val newVal =
            ifield.min(
              me.get(i,j)(x),
              ifield.+(me.get(i,k)(x), me.get(k,j)(x)))
          x = me.update(i, j, newVal)(x)
        }

    for (i <- 0 to 2 * nOfVars - 1)
      for (j <- 0 to 2 * nOfVars - 1) {
        val newVal =
          ifield.min(
            me.get(i, j)(x),
            ifield.+(
              me.get(i, signed(i))(x),
              me.get(signed(j), j)(x)))
      }
    x
  }

}
