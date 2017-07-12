package it.unich.jandom.domains.numerical.octagon

import scala.language.higherKinds
import scalaz._

trait Matrix[M[_]] // extends Functor[M] with Foldable[M] {
  extends Foldable[M] {

  def update[A](i: Int, j: Int, x: A)(m: M[A]): M[A]
  def update[A](updater: (Int, Int) => A)(m: M[A]): M[A]
  def get[A](i: Int, j: Int)(m: M[A]): Option[A]
  def combine[A, B, C](f: (A, B) => C)(ma: M[A], mb: M[B]): M[C]
}
