package it.unich.jandom.domains.numerical.octagon

import scala.language.higherKinds

/**
  * Created by fsestini on 7/11/17.
  *
  * Type of abstract octagons (elements of the octagon abstract domain),
  * parameterized by the underlying DBM and the type of elements.
  *
  * Extremely temporary...
  */
case class AbstractOctagon[M[_], A](dbm: M[A], e: DifferenceBoundMatrix[M]) {
  def join(other: AbstractOctagon[M, A])
          (implicit ifield: InfField[A], c: e.LatticeConstraint[A]): AbstractOctagon[M, A] =
    new AbstractOctagon(e.union(e.strongClosure(dbm), e.strongClosure(other.dbm)), e) // check if troo

  def meet(other: AbstractOctagon[M, A])
          (implicit c: e.LatticeConstraint[A]) : AbstractOctagon[M, A] =
    new AbstractOctagon(e.intersection(dbm, other.dbm), e) // check if troo

  def forget(): AbstractOctagon[M, A] = ???
  def assignment(): AbstractOctagon[M, A] = ???
  // etc...

  def top(implicit ifield: InfField[A]) =
    new AbstractOctagon(e.topDBM[A], e: DifferenceBoundMatrix[M])
  def bottom =
    new AbstractOctagon(e.bottomDBM[A], e: DifferenceBoundMatrix[M])
}
