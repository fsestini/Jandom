package it.unich.jandom.domains.numerical.octagon

import scala.language.higherKinds
import scalaz._

trait Matrix[M[_]] extends Functor[M] with Foldable[M] {
  def update[A](i: Int, j: Int, x: A)(m: M[A]): M[A]
  def update[A](updater: (Int, Int) => Option[A])(m: M[A]): M[A]
  def get[A](i: Int, j: Int)(m: M[A]): Option[A]
  // def empty[A](dimension: Int): M[A]

  // These should be in some typeclass instance
  def combine[A, B, C](f: (A, B) => C)(ma: M[A], mb: M[B]): M[C]
  // def toList[A](m: M[A]): List[A]
}
//
//object Matrix {
//  def fromFun[A, M[_]](dimension: Int, f: (Int, Int) => Option[A])(implicit e: Matrix[M]): M[A] = {
//    e.update(f)(e.empty(dimension))
//  }
//}