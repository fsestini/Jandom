package it.unich.jandom.domains.numerical.octagon.fast

import it.unich.jandom.domains.numerical._
import it.unich.jandom.domains.numerical.octagon._
import it.unich.jandom.domains.numerical.octagon.variables.CountOps._
import variables._
import VarIndexOps._
import org.scalatest.PropSpec
import org.scalatest.prop.PropertyChecks
import org.scalacheck._
import org.scalacheck.Arbitrary.arbitrary
import spire.math.Rational
import spire.math.RationalAlgebra
import it.unich.jandom.utils.numberext.RationalExt
import scala.language.higherKinds
import scala.reflect.ClassTag
import scala.language.reflectiveCalls

trait AbstractFastSpecification extends AbstractOctagonSpecification {

  type MainM[_]
  type SubM[_]

  type M[S,X] = CFastDBM[MainM, SubM, S, X]

  val mev: MEvidence[MainM, SubM]

  // Given a list of variables, set all entries involving variables in the list
  // to a non-infinite value, and leave all the others unchanged. Useful since
  // m' = linkVars(linkVars(m, xs), ys) where xs and ys are disjoint is such
  // that computeComponents(m') = { xs, ys }
  def linkVars(m: MainM[A], ixs: List[VarIndex])
    (implicit ifield: InfField[A]): MainM[A] = {

    def concerns(i: Int)(v: VarIndex): Boolean =
      varPlus(v) == i || varMinus(v) == i

    mev.ds.update((i, j) =>
      if (ixs.exists(concerns(i)) && ixs.exists(concerns(j))) ifield.zero
      else mev.ds.get(i, j)(m))(m)
  }

  def infinity(vc: VarCount)(implicit ifield: InfField[A]): MainM[A] =
    mev.dec.pure(vc, ifield.infinity)

  def GenSubsetOf(min: Int, max: Int): Gen[List[Int]] =
    if (min > max) Gen.const(Nil) else for {
      b <- Gen.oneOf(0, 1)
      xs <- GenSubsetOf(min + 1, max)
    } yield if (b == 1) min :: xs else xs

  def GenPartitionOf(list: List[Int]): Gen[List[List[Int]]] = list match {
    case Nil => Gen.const(Nil :: Nil)
    case x :: Nil => Gen.const((x :: Nil) :: Nil)
    case x :: xs => for {
      pss <- GenPartitionOf(xs)
      b <- Gen.oneOf(0,1)
    } yield pss match {
      case (p :: ps) => if (b == 1) ((x :: p) :: ps) else ((x :: Nil) :: p :: ps)
    }
  }

  property ("Decomposable matrices should get decomposed after closure") {
    forAll(Gen.choose(5,10)) { (n: Int) =>
      val vcount = n // 5 + n
      forAll(GenSubsetOf(0, vcount)) { (sub: List[Int]) =>
        forAll(GenPartitionOf(sub)) { (part: List[List[Int]]) =>

          // hack apparentemente necessario, visto che aggiungere suchThat(...)
          // ai generatori non ha alcun effetto
          val ppartitions = part.map(_.map(VarIndex))
          val partitions =
            if (ppartitions.length >= 2) ppartitions
            else allVars(VarCount(vcount)).toList :: Nil

          val inf = infinity(VarCount(vcount))
          val linked = partitions.foldLeft(inf)((m, p) => linkVars(m, p))
          val ncFast = FullDBM(linked, mev)

          val cFast: M[Closed, A] = ncFast.strongClosure(mev, ifield)
          cFast match {
            case CFast(DecomposedDBM(_, comps, _)) =>
              val set1: Set[Set[VarIndex]] = partitions.map(_.toSet).toSet
              val set2: Set[Set[VarIndex]] = comps.map(_.toSet).toSet
              set1 == set2
            case _ => false
          }
        }
      }
    }
  }


}
