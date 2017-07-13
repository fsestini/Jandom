package it.unich.jandom.domains.numerical.octagon

import scala.language.higherKinds

/**
  * Created by fsestini on 7/10/17.
  */

sealed trait Ordering
object GT extends Ordering
object LT extends Ordering
object EQ extends Ordering

trait Poset [P] {
  def compare(x: P, y: P): Option[Ordering]
}

trait Poset1[P[_]] {
  type PosetConstraint[_]
  def compare[A](x: P[A], y: P[A])(implicit evidence: PosetConstraint[A]): Option[Ordering]
}

trait Lattice1[L[_]] extends Poset1[L] {
  type LatticeConstraint[_]
  def union[A](x: L[A], y: L[A])(implicit evidence: LatticeConstraint[A]): L[A]
  def intersection[A](x: L[A], y: L[A])(implicit evidence: LatticeConstraint[A]): L[A]
}

trait CompleteLattice1[L[_]] extends Lattice1[L] {
  def top[A]: L[A]
  def bottom[A]: L[A]
}
