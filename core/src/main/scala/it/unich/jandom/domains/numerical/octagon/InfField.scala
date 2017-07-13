package it.unich.jandom.domains.numerical.octagon

import breeze.linalg.norm
import breeze.math.{Field, Ring}

/**
  * Created by fsestini on 7/10/17.
  */
trait InfField[@specialized(Double) A] extends Field[A] with Poset[A] {
  def infinity: A
  def max(x: A, y: A): A
  def min(x: A, y: A): A
  def half(x: A): A
}

object InfField {
  implicit object infFieldDouble extends InfField[Double] {
    def infinity: Double = Double.PositiveInfinity

    def -(a: Double, b: Double): Double = a - b

    def %(a: Double, b: Double): Double = a % b

    implicit val normImpl: norm.Impl[Double, Double] = new norm.Impl[Double, Double] {
      override def apply(v: Double): Double = v.abs
    }

    def zero: Double = 0
    def one: Double = 1

    def +(a: Double, b: Double): Double = a + b
    def *(a: Double, b: Double): Double = a * b
    def ==(a: Double, b: Double): Boolean = a == b
    def !=(a: Double, b: Double): Boolean = a != b
    def /(a: Double, b: Double): Double = a / b
    def pow(a: Double, b: Double): Double = pow(a, b)

    override def max(x: Double, y: Double): Double = math.max(x, y)
    override def min(x: Double, y: Double): Double = math.min(x, y)
    override def half(x: Double): Double = x / 2

    def compare(x: Double, y: Double): Option[Ordering] =
      (x < y, x == y) match {
        case (true, _) => Some(LT)
        case (_, true) => Some(EQ)
        case _ => Some(GT)
    }
  }
}
