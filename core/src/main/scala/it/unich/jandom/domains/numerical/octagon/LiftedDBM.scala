package it.unich.jandom.domains.numerical.octagon

import scalaz._
import std.option._, std.list._

/**
  * Created by fsestini on 7/10/17.
  */

sealed trait LiftedDBM[M <: DifferenceBoundMatrix[M]] extends DifferenceBoundMatrix[LiftedDBM[M]]

case class BottomDBM[M](nOfVariables: Int) extends LiftedDBM[M] {
  override def strongClosure(): LiftedDBM[M] = BottomDBM(nOfVariables)

  def incrementalClosure(): LiftedDBM[M] = ???

  def union(that: LiftedDBM[M]): LiftedDBM[M] = that

  def intersection(that: LiftedDBM[M]): LiftedDBM[M] = this

  def tryCompareTo[B >: LiftedDBM[M]](that: B)
    (implicit evidence$1: (B) => PartiallyOrdered[B]): Option[Int] =
      that match {
      case other: BottomDBM[M] => Some(0)
      case other: ActualDBM[M] => Some(-1)
      case _ => None
    }

}

case class ActualDBM[M <: DifferenceBoundMatrix[M]](m: M, nOfVariables: Int)
    extends LiftedDBM[M] {

  def any[A](l: List[A], p: A => Boolean): Boolean =
    l.foldRight(false)((x, y) => p(x) || y)

  def strongClosure(): LiftedDBM[M] = {
    val closed: M = m.strongClosure()
    val l: Option[List[Double]] =
      Applicative[Option].sequence(
        for (i <- 2 * nOfVariables) yield e.get(i, i)(closed))
    l match {
      case Some(list) => any(list, (x: Double) => x < 0) match {
        case true => new BottomDBM[M](nOfVariables)
        case false => {
          val fun: (Int, Int) => Option[Double] = (i, j) => i == j match {
            case true => Some(0)
            case false => m.matrix // e.get(i,j)(m)
          }
          new ActualDBM[M](e.update(fun)(m), nOfVariables)
        }
      }
      case None => throw new IllegalArgumentException()
    }
  }

  def incrementalClosure(): LiftedDBM[M] = ???

  def union(that: LiftedDBM[M]): LiftedDBM[M] = ???

  def intersection(that: LiftedDBM[M]): LiftedDBM[M] = ???

  def tryCompareTo[B >: LiftedDBM[M]](that: B)(implicit evidence$1: (B) => PartiallyOrdered[B]): Option[Int] = ???

}
//
//object LiftedMatrixInstance {
//  def LiftedMatrixIsMatrix[M <: DifferenceBoundMatrix[M]]: Matrix[({
//    type L[a] = LiftedDBM[M]})#L] = new Matrix[({type L[a] = LiftedDBM[M]})#L] {
//
//    def update[A](i: Int, j: Int, x: A)(m: LiftedDBM[M]): LiftedDBM[M] = m match {
//      case BottomDBM(n) => BottomDBM(n)
//      case ActualDBM(m, n) => ActualDBM(m.matrixEvidence.update(i, j, x)(m), n)
//    }
//
//    def update[A](updater: (Int, Int) => Option[A])(m: LiftedDBM[M]): LiftedDBM[M] = m match {
//      case BottomDBM(n) => BottomDBM(n)
//      case ActualDBM(m, n) => ActualDBM(m.matrixEvidence.update(updater)(m), n)
//    }
//
//    def get[A](i: Int, j: Int)(m: LiftedDBM[M]): Option[A] = m match {
//      case BottomDBM(n) => None
//      case ActualDBM(m, n) => m.matrixEvidence.get(i, j)(m)
//    }
//
//    def combine[A, B, C](f: (A, B) => C)(ma: LiftedDBM[M], mb: LiftedDBM[M]): LiftedDBM[M] =
//      ma match {
//        case BottomDBM(n) => BottomDBM(n)
//        case ActualDBM(m1, n1) => mb match {
//          case BottomDBM(n2) => BottomDBM(n2)
//          case ActualDBM(m2, _) => ActualDBM(m1.matrixEvidence.combine(f)(m1, m2), n1)
//        }
//      }
//
//    def toList[A](m: LiftedDBM[M]): List[A] = m match {
//      case BottomDBM(_) => Nil
//      case ActualDBM(m, _) => m.matrixEvidence.toList(m)
//    }
//  }
//}
//
