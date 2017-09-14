package it.unich.jandom.domains.numerical.octagon

import scala.language.higherKinds
import scalaz._
import variables._

trait Matrix[M[_]] // extends Functor[M] with Foldable[M] {
  extends Foldable[M] {

  def update[A](i: Int, j: Int, x: A)(m: M[A]): M[A]
  def update[A](updater: (Int, Int) => A)(m: M[A]): M[A]
  def get[A](i: Int, j: Int)(m: M[A]): A
  def combine[A, B, C](f: (A, B) => C)(ma: M[A], mb: M[B]): M[C]
  def pure[A](dimension: Dimension, x: A): M[A]
  def make[A](f: ((Int, Int) => A), dimension: Dimension): M[A]
  def dimension[A](m: M[A]): Dimension
}
