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
