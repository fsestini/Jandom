package it.unich.jandom.domains.numerical.octagon

sealed trait Ordering
case class GT() extends Ordering
case class LT() extends Ordering
case class EQ() extends Ordering

trait Poset1[P[_]] {
  def compare[A](x: P[A], y: P[A]): Option[Ordering]
}

trait Lattice1[L[_]] extends Poset1[L] {
  def union[A](x: L[A], y: L[A]): L[A]
  def intersection[A](x: L[A], y: L[A]): L[A]
}

trait DifferenceBoundMatrix[M[_]] extends Matrix[M] with Lattice1[M] {
    // extends Lattice[M] with PartiallyOrdered[M] {
  def strongClosure[A](m: M[A]): M[A]
  def incrementalClosure[A](m: M[A]): M[A]
  // def nOfVariables[A](m: M[A]): Int
}

