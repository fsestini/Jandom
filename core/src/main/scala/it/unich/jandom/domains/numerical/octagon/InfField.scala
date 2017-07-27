package it.unich.jandom.domains.numerical.octagon

import breeze.linalg.norm
import breeze.math.{Field, Ring}
import it.unich.jandom.utils.numberext.RationalExt

/**
  * Created by fsestini on 7/10/17.
  */
trait InfField[@specialized(Double) A] extends Field[A] with Toset[A] {
  def infinity: A
  def max(x: A, y: A): A
  def min(x: A, y: A): A
  def half(x: A): A
  def double(x: A): A
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
    def pow(a: Double, b: Double): Double = math.pow(a, b)

    override def max(x: Double, y: Double): Double = math.max(x, y)
    override def min(x: Double, y: Double): Double = math.min(x, y)
    override def half(x: Double): Double = x / 2
    override def double(x: Double): Double = x * 2

    def compare(x: Double, y: Double): Ordering =
      (x < y, x == y) match {
        case (true, _) => LT
        case (_, true) => EQ
        case _ => GT
    }
  }

  implicit object ifieldRationalExt extends InfField[RationalExt] {
    def infinity: RationalExt = RationalExt.PositiveInfinity
    def max(x: RationalExt, y: RationalExt): RationalExt = x.max(y)
    def min(x: RationalExt, y: RationalExt): RationalExt = x.min(y)
    def half(x: RationalExt): RationalExt = x / 2
    def double(x: RationalExt): RationalExt = x * 2
    def /(a: RationalExt, b: RationalExt): RationalExt = a / b
    def pow(a: RationalExt, b: RationalExt): RationalExt = ???
    def -(a: RationalExt, b: RationalExt): RationalExt = a - b
    def %(a: RationalExt, b: RationalExt): RationalExt = ???
    implicit val normImpl: norm.Impl[RationalExt, Double] = new norm.Impl[RationalExt, Double] {
      def apply(v: RationalExt): Double = v.abs.toDouble
    }
    def compare(x: RationalExt, y: RationalExt): Ordering =
      if (x == y) EQ else if (x < y) LT else GT
    def zero: RationalExt = 0
    def one: RationalExt = 1
    def +(a: RationalExt, b: RationalExt): RationalExt = a + b
    def *(a: RationalExt, b: RationalExt): RationalExt = a * b
    def ==(a: RationalExt, b: RationalExt): Boolean = a == b
    def !=(a: RationalExt, b: RationalExt): Boolean = a != b
  }
}
