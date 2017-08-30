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

  property ("T <> _|_ for different dim") {
    forAll (GenSmallInt) {
      (d: Int) => {
        forAll (GenSmallInt) {
          (delta: Int) => {
            val d1 = d + delta
            val topoct =  AbstractOctagon[OctagonDomain[FastDBM], FastDBM](e.topDBM[Double](VarCount(d)), oct, e)
            val botoct = AbstractOctagon[OctagonDomain[FastDBM], FastDBM](e.bottomDBM[Double](VarCount(d1)), oct, e)
            (
              !(topoct >= botoct) &
                !(botoct >= topoct)
            )
          }
        }
      }
    }
  }

  property ("T <> T for different dim") {
    forAll (GenSmallInt) {
      (d: Int) => {
        forAll (GenSmallInt) {
          (delta: Int) => {
            val d1 = d + delta
            val topoct =  AbstractOctagon[OctagonDomain[FastDBM], FastDBM](e.topDBM[Double](VarCount(d)), oct, e)
            val topoct2 = AbstractOctagon[OctagonDomain[FastDBM], FastDBM](e.topDBM[Double](VarCount(d1)), oct, e)
            (
              !(topoct >= topoct2) &
                !(topoct2 >= topoct)
            )
          }
        }
      }
    }
  }

  property ("_|_ <> _|_ for different dim") {
    forAll (GenSmallInt) {
      (d: Int) => {
        forAll (GenSmallInt) {
          (delta: Int) => {
            val d1 = d + delta
            val botoct =  AbstractOctagon[OctagonDomain[FastDBM], FastDBM](e.bottomDBM[Double](VarCount(d)), oct, e)
            val botoct2 = AbstractOctagon[OctagonDomain[FastDBM], FastDBM](e.bottomDBM[Double](VarCount(d1)), oct, e)
            (!(botoct >= botoct2) &
              (botoct2 >= botoct)
            )
          }
        }
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

}
