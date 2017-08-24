package it.unich.jandom.domains.numerical.octagon

import VarIndexOps._
import CountOps._

// FunMatrix-based raw DBM implementation
sealed trait FunDBM[S, A] {
  def liftFromInner(f: FunMatrix[A] => FunMatrix[A])(implicit ifield: InfField[A]): FunDBM[S, A]
  def union(other: FunDBM[S, A])(implicit infField: InfField[A]): FunDBM[S, A]
  def decideState: DBMIxed[FunDBM, A]
  val innerMatrix: Option[FunMatrix[A]]
  def noOfVariables: VarCount
}
object FunDBMInstance {
  implicit val funDBM: DifferenceBoundMatrix[FunDBM] { type PosetConstraint[A] = InfField[A] } = ???
}

// This is the simplified strong closure algorithm from
// Bagnara et al., Widening Operators for Weakly-Relational Numeric Abstractions.
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

  def nullCheck[A](m: FunMatrix[A])(implicit ifield: InfField[A]): Option[FunMatrix[A]] = {
    val negative: Boolean = allIndices(m.dimension).exists((i) =>
      ifield.compare(me.get(i, i)(m), ifield.zero) == LT)
    if (negative) None else {
      val updater: (Int, Int) => A = (i, j) =>
        if (i == j) ifield.zero else me.get(i, j)(m)
      Some(me.update(updater)(m))
    }
  }

  def strengthen[A](dbm: FunMatrix[A])(implicit ifield: InfField[A]): FunMatrix[A] =
    grid(dbm.dimension)
      .foldLeft(dbm)((x, pair) => pair match {
        case (i, j) => {
          val newVal =
            ifield.min(
              me.get(i, j)(x),
              ifield.half(ifield.+(me.get(i, signed(i))(x), me.get(signed(j), j)(x))))
          me.update(i, j, newVal)(x)
        }
      })

  def strongClosure[A](dbm: FunMatrix[A])(implicit ifield: InfField[A]): Option[FunMatrix[A]] = {
    val closed: FunMatrix[A] = (for {
      k <- allIndices(dbm.dimension)
      (i, j) <- grid(dbm.dimension)
    } yield (k, i, j)).foldLeft(dbm)((x, triple) => triple match {
      case (k, i, j) => {
        val newVal =
          ifield.min(
            me.get(i,j)(x),
            ifield.+(me.get(i,k)(x), me.get(k,j)(x)))
        me.update(i, j, newVal)(x)
      }
    })

    nullCheck(strengthen(closed))
  }

  def incrementalClosure[A](vi: VarIndex)(dbm: FunMatrix[A])
                           (implicit ifield: InfField[A]): Option[FunMatrix[A]] = {
    val p: ((Int, Int, Int)) => Boolean = {
      case (k, i, j) =>
        k == vi.i || k == signed(vi.i) ||
          i == vi.i || i == signed(vi.i) ||
          j == vi.i || j == signed(vi.i)
    }

    val iclosed: FunMatrix[A] =
      (for { k <- allIndices(dbm.dimension) ;
             (i, j) <- grid(dbm.dimension) } yield (k, i, j))
        .filter(p)
        .foldLeft(dbm)((x, triple) => triple match {
          case (k, i, j) => {
            val newVal =
              ifield.min(me.get(i,j)(x), ifield.+(me.get(i,k)(x), me.get(k,j)(x)))
            me.update(i, j, newVal)(x)
          }
        })

    nullCheck(strengthen(iclosed))
  }

}

object VarMapping {

  def varMapImageSize(f: VarIndex => Option[VarIndex], nOfVars: VarCount): VarCount =
    VarCount(
      allVars(nOfVars).map((v) => f(v) match {
        case Some(_) => 1
        case None => 0
      }).sum)

  def mapVariablesAux[A](f: (VarIndex) => Option[VarIndex], nOfVars: VarCount)
                        (m: FunMatrix[A]): FunMatrix[A] = {

    def inverse(f: VarIndex => Option[VarIndex])(v: VarIndex): VarIndex = {
      val vars: Seq[VarIndex] = allVars(nOfVars)
      val opts: Seq[Option[VarIndex]] = vars.map(x => f(x) match {
        case Some(vres) => if (v == vres) Some(x) else None
        case None => None
      })
      val opt: Option[VarIndex] = opts.flatten.headOption
      opt match {
        case Some(vres) => vres
        case None => throw new IllegalArgumentException("inverse out of bounds")
      }
    }

    def inverseIx(f: VarIndex => Option[VarIndex])(i: Int): Int =
      toIndexAndCoeff(i) match {
        case (v, Positive) => varPlus(inverse(f)(v))
        case (v, Negative) => varMinus(inverse(f)(v))
      }

    val newCount: VarCount = varMapImageSize(f, nOfVars)
    FunMatrix((i, j) => {
      m(inverseIx(f)(i), inverseIx(f)(j))
    }, doubledVarCount(newCount))
  }
}
