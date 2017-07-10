package it.unich.jandom.domains.numerical.octagon

import scala.language.higherKinds
import scalaz._
import std.option._, std.list._

class DenseDBM[M[_], A](m: M[A], nOfVariables: Int, e: Matrix[M])
    extends DifferenceBoundMatrix[DenseDBM[M, _]] {
  // val nni: Int

  private def cloTransform(k: Int)(m: M[Double]): M[Double] = {
    val f: (Int, Int) => Option[Double] = (i: Int, j: Int) => {
      val l: List[Option[Double]] = List(
        e.get(i,j)(m),
        Apply[Option].apply2(e.get(i, 2 * k - 1)(m), e.get(2 * k - 1, j)(m))(_ + _),
        Apply[Option].apply2(e.get(i, 2 * k)(m), e.get(2 * k, j)(m))(_ + _),
        Apply[Option].apply2(
          e.get(i, 2 * k - 1)(m),
          Apply[Option].apply2(
            e.get(2 * k - 1, 2 * k)(m),
            e.get(2 * k, j)(m))(_ + _)
        )(_ + _),
        Apply[Option].apply2(
          e.get(i, 2 * k)(m),
          Apply[Option].apply2(
            e.get(2 * k, 2 * k - 1)(m),
            e.get(2 * k - 1, j)(m))(_ + _)
        )(_ + _)
      )
      Applicative[Option].sequence(l) match {
        case Some(resList) => Some(resList.fold(Double.PositiveInfinity)(math.min _))
        case None => None
      }
    }
    e.update(f)(m)
  }

  private def close(m: M[Double]): M[Double] = {
    var x: M[Double] = m
    for (k <- 1 to nOfVariables) {
      x = cloTransform(k)(x)
    }
    x
  }

  private def signed(i: Int): Int = if (i % 2 != 0) (i + 1) else (i - 1)

  private def strengthen(m: M[Double]): M[Double] = {
    val updater: (Int, Int) => Option[Double] = (i: Int, j: Int) => for {
          a <- e.get(i, j)(m)
          b <- e.get(i, signed(i))(m)
          c <- e.get(signed(j), j)(m)
        } yield math.min(a, (b + c) / 2)
    e.update(updater)(m)
  }

  def strongClosure() = new DenseDBM(strengthen(close(m)), nOfVariables, e)

  def incrementalClosure() = ???

  def union(that: DenseDBM[M]) =
    new DenseDBM[M](
      e.combine((x: Double, y: Double) => math.max(x,y))(m, that.m),
      nOfVariables, e)

  def intersection(that: DenseDBM[M]): DenseDBM[M] =
    new DenseDBM(
      e.combine((x: Double, y: Double) => math.min(x,y))(m, that.m),
      nOfVariables, e)

  def all[A](l: List[A], p: A => Boolean): Boolean =
    l.foldRight(true)((x, y) => p(x) && y)

  def tryCompareTo[B >: DenseDBM[M]](that: B)
           (implicit evidence$1: B => PartiallyOrdered[B]): Option[Int] =
    that match {
      case other: DenseDBM[M] => nOfVariables == other.nOfVariables match {
        case true => {
          val qqq: M[Option[Int]] = e.combine(
            (x: Double, y: Double) => tryCompareTo(x, y))(m, other.m)
          val mm: List[Option[Int]] = e.toList(qqq)
          Applicative[Option].sequence(mm) match {
            case Some(list) => {
              (all(list, _ == -1), all(list, _ == 0), all(list, _ == 1)) match {
                case (true, _, _) => Some(-1)
                case (_, true, _) => Some(0)
                case (_, _, true) => Some(1)
                case _ => None
              }
            }
            case None => None
          }
        }
        case false => None
      }
      case _ => None
    }

  def matrixEvidence: Matrix[({type L[a] = DenseDBM[M]})#L] = new Matrix[({type L[a] = DenseDBM[M]})#L] {
    def update[A](i: Int, j: Int, x: A)(m: DenseDBM[M]): DenseDBM[M] =
      new DenseDBM(e.update(i, j, x)(m.m), nOfVariables, e)

    def update[A](updater: (Int, Int) => Option[A])(m: DenseDBM[M]): DenseDBM[M] = ???

    def get[A](i: Int, j: Int)(m: DenseDBM[M]): Option[A] = ???

    def combine[A, B, C](f: (A, B) => C)(ma: DenseDBM[M], mb: DenseDBM[M]): DenseDBM[M] = ???

    def toList[A](m: DenseDBM[M]): List[A] = ???
  }
}

object DenseDBM {
  def apply(nOfVariables: Int)(implicit evidence: Matrix[FunMatrix]): DenseDBM[FunMatrix] =
    new DenseDBM[FunMatrix](
      new FunMatrix[Double](
        (_,_) => Some(Double.PositiveInfinity),
        2 * nOfVariables),
      nOfVariables, evidence)
}
