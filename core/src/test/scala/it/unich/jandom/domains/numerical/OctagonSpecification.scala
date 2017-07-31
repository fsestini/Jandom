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

package it.unich.jandom.domains.numerical.octagon
import it.unich.jandom.domains.numerical._
import org.scalatest.PropSpec
import org.scalatest.prop.{PropertyChecks, Checkers}
import org.scalacheck.Gen
import org.scalacheck.Arbitrary
import Arbitrary._

class OctagonSpecification extends PropSpec with PropertyChecks {
  val e = FunDBMInstance.funDBM
  val oct = new OctagonDomain(e)
  val box = BoxDoubleDomain(false)

  def GenSmallInt : Gen[Int] = Gen.choose(1, 16)
  def GenSmallEvenInt : Gen[Int] = for (n <- GenSmallInt) yield (n * 2)
  def GenInf : Gen[Double] = Gen.const(Double.PositiveInfinity)
  def GenNonnegDoubles : Gen[Double] = Gen.choose(0, Double.MaxValue)
  def GenDoublesAndInf(p: Int = 10) : Gen[Double] = Gen.frequency(
    (1, GenInf),
    (p-1, GenNonnegDoubles)
  )

  def GenOrderedPair : Gen[(Double, Double)] = for {
    low <- arbitrary[Double]
    high <- Gen.choose(low, Double.MaxValue)
  } yield (low, high)

  def GenOrderedDistinctPair : Gen[(Double, Double)] = for {
    low <- arbitrary[Double]
    high <- Gen.choose(low, Double.MaxValue).suchThat(_ > low)
  } yield (low, high)

  def GenMatrix : Gen[FunMatrix[Double]] = for {
    d <- GenSmallEvenInt
    rowSeq <- Gen.containerOfN[Array, Double](d, GenDoublesAndInf())
    arrayOfRows <- Gen.containerOfN[Array, Array[Double]](d, rowSeq)
  } yield (
      new FunMatrix[Double](
        (i: Int, j: Int) => arrayOfRows(i)(j), d))

  implicit def arbBox : Arbitrary[box.Property] =
    Arbitrary {

      for {
        n <- Gen.choose(1,20)
        pairs : Array[(Double, Double)] <- Gen.containerOfN[Array, (Double, Double)](n, GenOrderedPair)
      } yield (new box.Property(pairs.unzip._1, pairs.unzip._2, false))
    }

  property("Can create from interval") {
    forAll {
      (b1: box.Property)
      =>
      val a = AbstractOctagon.fromInterval(b1, oct, e)
      true
    }
  }

  // We are *probably* no longer interested in this?
  //
  // property ("After trivial singlePositiveExactAssignment to C toInterval is [C,C]") = forAll {
  //   (c: Int) => {
  //     val a = AbstractOctagon(FunDBMInstance.funDBM.topDBM[Double](1), oct, e)
  //     val b = AbstractOctagon(DBMUtils.singlePositiveExactAssignment(VarIndex(0), c.toDouble)(a.dbm, FunDBMInstance.funDBM), oct, e)
  //     b.toInterval.high.head == c &&
  //     b.toInterval.low.head == c
  //   }
  // }

  property ("T.{x <- C} toInterval is [C,C]") {
    forAll {
      (c: Int) => {
        val a = AbstractOctagon(FunDBMInstance.funDBM.topDBM[Double](1), oct, e)
        val b = a.linearAssignment(0, LinearForm.c(c))
        (b.toInterval.high.head == c &
         b.toInterval.low.head == c)
      }
    }
  }

  property("toInterval.fromInterval == id") {
    forAll {
      (b1: box.Property)
      =>
      AbstractOctagon.fromInterval(b1, oct, e).toInterval == b1
    }
  }

  property ("T >= _|_") {
    forAll (GenSmallInt) {
      case (d: Int) => {
        val topoct =  AbstractOctagon(e.topDBM[Double](d), oct, e)
        val botoct = AbstractOctagon(e.bottomDBM[Double](d), oct, e)
        (topoct >= botoct)
      }
    }
  }

  property ("T <> _|_ for different dim") {
    forAll (GenSmallInt) {
      case (d: Int) => {
        forAll (GenSmallInt) {
          case (delta: Int) => {
            val d1 = d + delta
            val topoct =  AbstractOctagon(e.topDBM[Double](d), oct, e)
            val botoct = AbstractOctagon(e.bottomDBM[Double](d1), oct, e)
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
      case (d: Int) => {
        forAll (GenSmallInt) {
          case (delta: Int) => {
            val d1 = d + delta
            val topoct =  AbstractOctagon(e.topDBM[Double](d), oct, e)
            val topoct2 = AbstractOctagon(e.topDBM[Double](d1), oct, e)
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
      case (d: Int) => {
        forAll (GenSmallInt) {
          case (delta: Int) => {
            val d1 = d + delta
            val botoct =  AbstractOctagon(e.bottomDBM[Double](d), oct, e)
            val botoct2 = AbstractOctagon(e.bottomDBM[Double](d1), oct, e)
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
      case (d: Int) => {
        val topoct =  AbstractOctagon(e.topDBM[Double](d), oct, e)
        val anotherTopoct = AbstractOctagon(e.bottomDBM[Double](d), oct, e)
        (topoct >= anotherTopoct)
      }
    }
  }

  property ("T U _|_ = T != _|_") {
    forAll (GenSmallInt) {
      case (d: Int) => {
        val topoct = AbstractOctagon(e.topDBM[Double](d), oct, e)
        val botoct = AbstractOctagon(e.bottomDBM[Double](d), oct, e)
        ((topoct union botoct) == topoct &
         (topoct union botoct) != botoct)
      }
    }
  }

  property ("T and _|_ intersection T = _|_") {
    forAll (GenSmallInt) {
      case (d: Int) => {
        val topoct = AbstractOctagon(e.topDBM[Double](d), oct, e)
        val botoct = AbstractOctagon(e.bottomDBM[Double](d), oct, e)
        ((topoct intersection botoct) != topoct &
         (topoct intersection botoct) == botoct)
      }
    }
  }

  property ("Union of [C1,C1], [C2, C2] == [C1, C2] w/C1 < C2") {
    forAll(GenOrderedDistinctPair) {
      case (c1: Double, c2: Double) => {
        val a = AbstractOctagon(e.topDBM[Double](1), oct, e)
        val b1 = a.linearAssignment(0, c1)
        val b2 = a.linearAssignment(0, c2)
        val union = b1 union b2
        (union.toInterval.isEmpty == false &
         union.toInterval.high.head == c2 &
         union.toInterval.low.head == c1)
      }
    }
  }

  property ("Intersection of [C1,C1], [C2, C2] == _|_ w/C1 < C2") {
    forAll(GenOrderedDistinctPair) {
      case (c1: Double, c2: Double) => {
        val a = AbstractOctagon(e.topDBM[Double](1), oct, e)
        val b1 = a.linearAssignment(0, c1)
        val b2 = a.linearAssignment(0, c2)
        val union = b1 intersection b2
        (union.isBottom == true &
         union.toInterval.isEmpty == true)
      }
    }
  }

  property ("[C1,C2] <= [C3<C1, C4>C2]") {
    forAll(GenOrderedDistinctPair) {
      case (x1: Double, x2: Double) =>
        forAll(GenOrderedDistinctPair) {
          case (x3: Double, x4: Double) => {
            val a = AbstractOctagon(e.topDBM[Double](1), oct, e)
            val c1 = LinearForm.c(x1);
            val c2 = LinearForm.c(x2);
            val posC = (x4 - x3) // x3 - x4 is certainly positive
            assert(posC > 0)
            val c3 = LinearForm.c(x1 - posC);
            val c4 = LinearForm.c(x2 + posC);
            val b1 = a.linearAssignment(0, c1)
            val b2 = a.linearAssignment(0, c2)
            val b3 = a.linearAssignment(0, c3)
            val b4 = a.linearAssignment(0, c4)
            val i1 = b1 union b2 // i1 = [x1,x2]
            val i2 = b3 union b4 // i2 = [x3,x4]
            (i1 <= i2)
          }
        }
      }
  }

  // TODO
  // property ("strongClosure(m) is closed") = ???

}
