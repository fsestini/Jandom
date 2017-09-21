/**
  * Copyright 2017 Davide Dal Bianco, Filippo Sestini, Tobia Tesan
  *
  * This file is part of JANDOM: JVM-based Analyzer for Numerical DOMains
  * JANDOM is free software: you can redistribute it and/or modify
  * it under the terms of the GNU General Public License as published by
  * the Free Software Foundation, either version 3 of the License, or
  * (at your option) any later version.
  *
  * JANDOM is distributed in the hope that it will be useful,
  * but WITHOUT ANY WARRANTY; without even the implied warranty of a
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  * GNU General Public License for more details.
  *
  * You should have received a copy of the GNU General Public License
  * along with JANDOM.  If not, see <http://www.gnu.org/licenses/>.
  */

package it.unich.jandom.domains.numerical.octagon.fast
import it.unich.jandom.domains.numerical.octagon._
import it.unich.jandom.domains.numerical._
import variables._
import VarIndexOps._
import CountOps._
import InfField._
import org.scalatest.PropSpec
import org.scalatest.prop.PropertyChecks
import spire.math.Rational
import spire.math.RationalAlgebra
import org.scalacheck.Gen
import org.scalacheck.Arbitrary
import org.scalacheck.Arbitrary.arbitrary
import it.unich.jandom.utils.numberext.RationalExt

class FastOctagonSpecification extends PropSpec with PropertyChecks {
  val mev = MEvidence(
    HalfMatrixDenseSparseInstance.halfMatrixDecomposableInstance,
    HalfMatrixDenseSparseInstance.halfMatrixDenseSparseInstance,
    HalfMatrixDenseSparseInstance.halfMatrixDecomposableInstance)
  val inf = Double.PositiveInfinity
  val e = CFDBMInstance.instance(mev)
  type FastDBM[S, A] = CFastDBM[HalfMatrix, HalfSubMatrix, S, A]
  val box =  BoxRationalDomain()
  val oct = new OctagonDomain[FastDBM, RationalExt, BoxRationalDomain](e, box)
  val Utils = new it.unich.jandom.domains.numerical.octagon.testutils.Utils(box)
  def id (x: RationalExt): RationalExt = x
  import Utils._
  type DOM = OctagonDomain[FastDBM, RationalExt, BoxRationalDomain]

  implicit def arbBox : Arbitrary[box.Property] =
    Arbitrary {
      for {
        n <- Gen.choose(1,20)
        pairs : Array[(RationalExt, RationalExt)] <- Gen.containerOfN[Array, (RationalExt, RationalExt)](n, GenOrderedPair)
      } yield (new box.Property(pairs.unzip._1, pairs.unzip._2, false))
    }

  //////////////////////////////////////////////////////////////////////////////
  // Begin properties
  //////////////////////////////////////////////////////////////////////////////

  property ("T.{x <- C} toInterval is [C,C]") {
    forAll {
      (c: Int) => {
        val top = e.topDBM[Double](VarCount(1))
        val a = AbstractOctagon[DOM, FastDBM, RationalExt, BoxRationalDomain](e.topDBM[RationalExt](VarCount(1)), oct, box, e)
        val b = a.linearAssignment(0, LinearForm.c(c))
        (b.toInterval.high.head == c &
         b.toInterval.low.head == c)
      }
    }
  }

  property("Can create from interval") {
    forAll {
      (b1: box.Property)
      =>
      val a = AbstractOctagon.fromInterval[DOM, FastDBM, RationalExt, BoxRationalDomain] (b1, oct, box, e)
      true
    }
  }


  property("toInterval.fromInterval == id") {
    forAll {
      (b1: box.Property)
      =>
      AbstractOctagon.fromInterval[DOM, FastDBM, RationalExt, BoxRationalDomain](b1, oct, box, e).toInterval == b1
    }
  }

  property ("T >= _|_") {
    forAll (GenSmallInt) {
      (d: Int) => {
        val topoct =  AbstractOctagon[DOM, FastDBM, RationalExt, BoxRationalDomain](e.topDBM[RationalExt](VarCount(d)), oct, box, e)
        val botoct = AbstractOctagon[DOM, FastDBM, RationalExt, BoxRationalDomain](e.bottomDBM[RationalExt](VarCount(d)), oct, box, e)
        (topoct >= botoct)
      }
    }
  }

  property ("T >= T") {
    forAll (GenSmallInt) {
      (d: Int) => {
        val topoct =  AbstractOctagon[DOM, FastDBM, RationalExt, BoxRationalDomain](e.topDBM[RationalExt](VarCount(d)), oct, box, e)
        val anotherTopoct = AbstractOctagon[DOM, FastDBM, RationalExt, BoxRationalDomain](e.bottomDBM[RationalExt](VarCount(d)), oct, box, e)
        (topoct >= anotherTopoct)
      }
    }
  }

  property ("T U _|_ = T != _|_") {
    forAll (GenSmallInt) {
      (d: Int) => {
        val topoct = AbstractOctagon[DOM, FastDBM, RationalExt, BoxRationalDomain](e.topDBM[RationalExt](VarCount(d)), oct, box, e)
        val botoct = AbstractOctagon[DOM, FastDBM, RationalExt, BoxRationalDomain](e.bottomDBM[RationalExt](VarCount(d)), oct, box, e)
        ((topoct union botoct) == topoct &
         (topoct union botoct) != botoct)
      }
    }
  }

  property ("T and _|_ intersection T = _|_") {
    forAll (GenSmallInt) {
      (d: Int) => {
        val topoct = AbstractOctagon[DOM, FastDBM, RationalExt, BoxRationalDomain](e.topDBM[RationalExt](VarCount(d)), oct, box, e)
        val botoct = AbstractOctagon[DOM, FastDBM, RationalExt, BoxRationalDomain](e.bottomDBM[RationalExt](VarCount(d)), oct, box, e)
        ((topoct intersection botoct) != topoct &
         (topoct intersection botoct) == botoct)
      }
    }
  }

  property ("Union of [C1,C1], [C2, C2] == [C1, C2] w/C1 < C2") {
    forAll(GenOrderedDistinctFinitePair) {
      (c : (Rational, Rational)) => {
        val a = AbstractOctagon[DOM, FastDBM, RationalExt, BoxRationalDomain](e.topDBM[RationalExt](VarCount(1)), oct, box, e)
        val b1 = a.linearAssignment(0, c._1)
        val b2 = a.linearAssignment(0, c._2)
        val union = b1 union b2
        assert(c._1 < c._2)
        (union.toInterval.isEmpty == false &
         union.toInterval.high.head == c._2 &
         union.toInterval.low.head == c._1)
      }
    }
  }

  property ("Intersection of [C1,C1], [C2, C2] == _|_ w/C1 < C2") {
    forAll(GenOrderedDistinctFinitePair) {
      (c : (Rational, Rational)) => {
        val a = AbstractOctagon[DOM, FastDBM, RationalExt, BoxRationalDomain](e.topDBM[RationalExt](VarCount(1)), oct, box, e)
        val b1 = a.linearAssignment(0, c._1)
        val b2 = a.linearAssignment(0, c._2)
        val intersection = b1 intersection b2
        assert(c._1 < c._2)
        (intersection.isBottom == true &
         intersection.toInterval.isEmpty == true)
      }
    }
  }

  property ("[C1,C2] <= [C3<C1, C4>C2]") {
    forAll(GenOrderedDistinctFinitePair) {
      (c: (Rational, Rational)) =>
        forAll(GenOrderedDistinctFinitePair) {
          (d: (Rational, Rational)) => {
            val c1 = c._1
            val c2 = c._2
            assert(c._1 < c._2)
            assert(d._1 < d._2)
            val posC = (d._2 - d._1) // is certainly positive
            assert(posC > 0)
            val c3 = c1 - posC
            val c4 = c2 + posC
            assert(c3 < c1 & c1 < c4 & c4 > c2)
            val a = AbstractOctagon[DOM, FastDBM, RationalExt, BoxRationalDomain](e.topDBM[RationalExt](VarCount(1)), oct, box, e)
            val b1 = a.linearAssignment(0, LinearForm.c(c1))
            val b2 = a.linearAssignment(0, LinearForm.c(c2))
            val b3 = a.linearAssignment(0, LinearForm.c(c3))
            val b4 = a.linearAssignment(0, LinearForm.c(c4))
            val i1 = b1 union b2 // i1 = [c1,c2]
            val i2 = b3 union b4 // i2 = [c3,c4]
            (i1 <= i2)
          }
        }
    }
  }

  // Given a list of variables, set all entries involving variables in the list
  // to a non-infinite value, and leave all the others unchanged. Useful since
  // m' = linkVars(linkVars(m, xs), ys) where xs and ys are disjoint is such
  // that computeComponents(m') = { xs, ys }
  def linkVars[A](m: HalfMatrix[A], ixs: List[VarIndex])
    (implicit ifield: InfField[A]): HalfMatrix[A] = {

    def concerns(i: Int)(v: VarIndex): Boolean =
      varPlus(v) == i || varMinus(v) == i

    m.update((i, j) =>
      if (ixs.exists(concerns(i)) && ixs.exists(concerns(j)))
        ifield.zero
      else m(i, j))
  }

  def infinity[A](vc: VarCount)(implicit ifield: InfField[A]): HalfMatrix[A] =
    new HalfMatrix(varCountToDim(vc), ifield.infinity)

  property ("Decomposable matrices should get decomposed after closure") {
    forAll(Gen.choose(5,10)) {
      (n: Int) =>
      val vcount = n // 5 + n
      forAll(GenSubsetOf(0, vcount)) {
        (sub: List[Int]) =>
        forAll(GenPartitionOf(sub)) {
          (part: List[List[Int]]) =>

          // hack apparentemente necessario, visto che aggiungere suchThat(...)
          // ai generatori non ha alcun effetto
          val ppartitions = part.map(_.map(VarIndex))
          val partitions =
            if (ppartitions.length >= 2) ppartitions
            else allVars(VarCount(vcount)).toList :: Nil

          val hm: HalfMatrix[RationalExt] = infinity(VarCount(vcount))
          val linked = partitions.foldLeft(hm)((m, p) => linkVars(m, p))
          val ncFast = FullDBM(linked, mev)

          val cFast: CFastDBM[HalfMatrix, HalfSubMatrix, Closed, RationalExt] =
            ncFast.strongClosure(mev, ifieldRationalExt)
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

  property ("Check that strongClosure for HalfMatrix DBMs is coherent") {
    forAll(GenSmallEvenInt) {
      (dd: Int) =>
      val d = (dd.abs + 1) * 2 // hack apparentemente necessario due to shrinking
      forAll (GenHalfMatrix(d)) {
        (m : HalfMatrix[RationalExt]) =>
        mev.ds.strongClosure(m) match {
          case None => true
          case Some(closed) =>
            forAll (
              Gen.zip(
                Gen.choose(0,d - 1),
                Gen.choose(0,d - 1))) {
              case (i, j) =>
                val ibar = signed(i)
                val jbar = signed(j)
                closed(jbar, ibar) == closed(i, j)
            }
        }
      }
    }
  }

  property ("Check that strongClosure for HalfMatrix DBMs is closed") {
    forAll(GenSmallEvenInt) {
      (dd: Int) =>
      val d = (dd.abs + 1) * 2 // hack apparentemente necessario due to shrinking
      forAll (GenHalfMatrix(d)) {
        (m : HalfMatrix[RationalExt]) =>
        mev.ds.strongClosure(m) match {
          case None => assert(true)
          case Some(c) =>
            forAll (
              Gen.zip(
                Gen.choose(0,d - 1),
                Gen.choose(0,d - 1),
                Gen.choose(0,d - 1))) {
              case (i, j, k) =>
                if (i == j)
                  c(i,j) == 0
                else
                  c(i,j) <= c(i, k) + c(k, j)
            }
        }
      }
    }
  }

  property ("Check that strongClosure for HalfMatrix DBMs has null diagonal") {
    forAll(GenSmallEvenInt) {
      (dd: Int) =>
      val d = (dd.abs + 1) * 2 // hack apparentemente necessario due to shrinking
      forAll (GenHalfMatrix(d)) {
        (m : HalfMatrix[RationalExt]) =>
        mev.ds.strongClosure(m) match {
          case None => assert(true)
          case Some(closed) =>
            (0 until d - 1).forall(i => closed(i, i) == 0)
        }
      }
    }
  }

}
