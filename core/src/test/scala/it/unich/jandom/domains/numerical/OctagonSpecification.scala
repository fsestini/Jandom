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
import VarIndexOps._
import org.scalatest.PropSpec
import org.scalatest.prop.PropertyChecks
import org.scalacheck.Gen
import org.scalacheck.Arbitrary
import org.scalacheck.Arbitrary.arbitrary
import org.scalacheck.Gen
import spire.math.Rational
import spire.math.RationalAlgebra

class OctagonSpecification extends PropSpec with PropertyChecks {
  val r = new RationalAlgebra()
  val e = FunDBMInstance.funDBM
  val oct = new OctagonDomain(e)
  val box = BoxDoubleDomain(false)

  def GenSmallInt : Gen[Int] = Gen.choose(1, 5)
  def GenSmallEvenInt : Gen[Int] = for (n <- GenSmallInt) yield (n * 2)
  def GenInf : Gen[Double] = Gen.const(Double.PositiveInfinity)
  def GenDoublesAndInf(pInf: Int) : Gen[Double] = Gen.frequency(
    (pInf, GenInf),
    (100 - pInf, Gen.choose(Float.MinValue.toDouble, Float.MaxValue.toDouble))
    // We generate only single precision floats to avoid false positives due to 754 edge cases
  )

  def GenChain : Gen[Seq[Double]] = Gen.containerOf[List, Double](arbitrary[Double]).map(_.sorted)

  implicit def arbRational: Arbitrary[Rational] =
    Arbitrary {
      for {
        n <- arbitrary[Int]
        d <- arbitrary[Int]
      } yield(r.fromInt(n)) / r.fromInt(Math.max(1,Math.abs(d))) // max(1,d) is a hackish way to avoid n/0
    }

  def GenArbitraryLf(n: Int): Gen[LinearForm] = for
  {
    coeffs <- Gen.containerOfN[List, Rational](n, arbitrary[Rational])
  } yield new DenseLinearForm(coeffs)

  def GenExactLf(n: Int): Gen[LinearForm] = for
  {
    vi <- Gen.choose(0, n - 1)
    ci <- Gen.oneOf(+1, -1)
    vj <- Gen.choose(0, n - 1)
    cj <- Gen.oneOf(+1, 0, -1)
  } yield new DenseLinearForm(
    Array.fill[Rational](n)(0)
      .updated(vi, r.toRational(ci))
      .updated(vj, r.toRational(cj))
      .toSeq
  )

  def GenLf(n: Int, pExact: Int = 10) : Gen[LinearForm] = Gen.frequency(
    (100 - pExact, GenArbitraryLf(n)),
    (pExact, GenExactLf(n))
  )

  // These are for disabling shrinking
  // From https://gist.github.com/davidallsopp/f65d73fea8b5e5165fc3
  //
  // TODO Find a better way

  import org.scalacheck.Shrink
  implicit val noShrink: Shrink[Int] = Shrink.shrinkAny
  implicit val noShrink2: Shrink[(Double, Double)] = Shrink.shrinkAny
  implicit val noShrink3: Shrink[FunMatrix[Double]] = Shrink.shrinkAny

  def GenOrderedPair : Gen[(Double, Double)] = for {
    low <- Gen.choose(Float.MinValue.toDouble, Float.MaxValue.toDouble)
    high <- Gen.choose(low, Float.MaxValue.toDouble)
    // We generate only single precision floats to avoid false positives due to 754 edge cases
  } yield (low, high)

  val GenOrderedDistinctPair = GenOrderedPair.suchThat((pair:(Double, Double)) =>pair._2 > pair._1)

  def GenFunMatrix(d: Int, pTop: Int = 20, pInf: Int = 20) : Gen[FunMatrix[Double]] = for {
    rowSeq <- Gen.containerOfN[Array, Double](d, GenDoublesAndInf(pInf))
    arrayOfRows <- Gen.containerOfN[Array, Array[Double]](d, rowSeq)
  } yield new FunMatrix[Double](
      ((i: Int, j: Int) =>
        if (i == j) 0
        else arrayOfRows(i)(j)), Dimension(d))

  def GenTop(nOfVars: Int): Gen[FunDBM[Closed, Double]] = Gen.const(new ClosedFunDBM(
    // Caveat: top is not strongly closed from the defn, but we have it as an instance of FunDBM[Closed,_].
        FunMatrix[Double]((i: Int, j: Int) =>
        if (i == j) 0 // If (i,i) != 0 we have bottom
        else Double.PositiveInfinity, Dimension(nOfVars*2)))
  )

  def GenClosedFunDBMOrTop(nOfVars: Int, pBot: Int = 10, pTop: Int = 20, pInf: Int = 30) : Gen[FunDBM[Closed, Double]] =
    Gen.frequency(
      (100 - pTop, GenTop(nOfVars)),
      (pTop, GenClosedFunDBM(nOfVars, pInf))
    )


  def GenClosedFunDBM(nOfVars: Int, pInf: Int = 30) : Gen[FunDBM[Closed, Double]] =
      for { m <- GenFunMatrix(nOfVars * 2, pInf) }
      yield
      {
        // caveat: there is no guarantee re: the distribution of bottoms, should probably include a few pre-computed non-bottom ones?
        val closure = BagnaraStrongClosure.strongClosure(m)
        if (closure == None) {
          new BottomFunDBM(VarCount(nOfVars))
        } else {
          new ClosedFunDBM(closure.get)
        }
      }

  implicit def arbBox : Arbitrary[box.Property] =
    Arbitrary {
      for {
        n <- Gen.choose(1,20)
        pairs : Array[(Double, Double)] <- Gen.containerOfN[Array, (Double, Double)](n, GenOrderedPair)
      } yield (new box.Property(pairs.unzip._1, pairs.unzip._2, false))
    }


  def checkIsLegal(m : FunMatrix[Double]) : Boolean =
    (1 until m.dimension.dim).forall(
      (i)=>(1 until m.dimension.dim).forall(
        (j)=>
        !m(i,j).isNaN
          &&
        (m(i,j) >= Double.MinValue
          &
          m(i,j) <= Double.MaxValue)
          |
          m(i,j) == Double.PositiveInfinity
      )
    )


  //////////////////////////////////////////////////////////////////////////////
  // Begin properties
  //////////////////////////////////////////////////////////////////////////////

  property("Can create from interval") {
    forAll {
      (b1: box.Property)
      =>
      val a = AbstractOctagon.fromInterval(b1, oct, e)
      true
    }
  }

  property ("T.{x <- C} toInterval is [C,C]") {
    forAll {
      (c: Int) => {
        val a = AbstractOctagon(FunDBMInstance.funDBM.topDBM[Double](VarCount(1)), oct, e)
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
      (d: Int) => {
        val topoct =  AbstractOctagon(e.topDBM[Double](VarCount(d)), oct, e)
        val botoct = AbstractOctagon(e.bottomDBM[Double](VarCount(d)), oct, e)
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
            val topoct =  AbstractOctagon(e.topDBM[Double](VarCount(d)), oct, e)
            val botoct = AbstractOctagon(e.bottomDBM[Double](VarCount(d1)), oct, e)
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
            val topoct =  AbstractOctagon(e.topDBM[Double](VarCount(d)), oct, e)
            val topoct2 = AbstractOctagon(e.topDBM[Double](VarCount(d1)), oct, e)
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
            val botoct =  AbstractOctagon(e.bottomDBM[Double](VarCount(d)), oct, e)
            val botoct2 = AbstractOctagon(e.bottomDBM[Double](VarCount(d1)), oct, e)
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
        val topoct =  AbstractOctagon(e.topDBM[Double](VarCount(d)), oct, e)
        val anotherTopoct = AbstractOctagon(e.bottomDBM[Double](VarCount(d)), oct, e)
        (topoct >= anotherTopoct)
      }
    }
  }

  property ("T U _|_ = T != _|_") {
    forAll (GenSmallInt) {
      (d: Int) => {
        val topoct = AbstractOctagon(e.topDBM[Double](VarCount(d)), oct, e)
        val botoct = AbstractOctagon(e.bottomDBM[Double](VarCount(d)), oct, e)
        ((topoct union botoct) == topoct &
         (topoct union botoct) != botoct)
      }
    }
  }

  property ("T and _|_ intersection T = _|_") {
    forAll (GenSmallInt) {
      (d: Int) => {
        val topoct = AbstractOctagon(e.topDBM[Double](VarCount(d)), oct, e)
        val botoct = AbstractOctagon(e.bottomDBM[Double](VarCount(d)), oct, e)
        ((topoct intersection botoct) != topoct &
         (topoct intersection botoct) == botoct)
      }
    }
  }

  property ("Union of [C1,C1], [C2, C2] == [C1, C2] w/C1 < C2") {
    forAll(GenOrderedDistinctPair) {
      (c : (Double, Double)) => {
        val a = AbstractOctagon(e.topDBM[Double](VarCount(1)), oct, e)
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
    forAll(GenOrderedDistinctPair) {
      (c : (Double, Double)) => {
        val a = AbstractOctagon(e.topDBM[Double](VarCount(1)), oct, e)
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
    forAll(GenOrderedDistinctPair) {
      (c: (Double, Double)) =>
        forAll(GenOrderedDistinctPair) {
          (d: (Double, Double)) => {
            val c1 = c._1
            val c2 = c._2
            assert(c._1 < c._2)
            assert(d._1 < d._2)
            val posC = (d._2 - d._1) // is certainly positive
            assert(posC > 0)
            val c3 = c1 - posC
            val c4 = c2 + posC
            assert(c3 < c1 & c1 < c4 & c4 > c2)
            val a = AbstractOctagon(e.topDBM[Double](VarCount(1)), oct, e)
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

  property ("Check that strongClosure is coherent (condition 1 of 3 for strong closure)") {
    forAll(GenSmallEvenInt) {
      (d: Int) =>
      forAll (GenFunMatrix(d)) {
        (m : FunMatrix[Double]) =>
        BagnaraStrongClosure.strongClosure(m) match {
          case None => false
          case Some(c) =>
            forAll (
              Gen.zip(
                Gen.choose(0,c.dimension.dim - 1),
                Gen.choose(0,c.dimension.dim - 1))) {
              case (i, j) =>
                val ibar = signed(i)
                val jbar = signed(j)
                c(jbar, ibar) == c(i, j)
            }
        }
      }
    }
  }

  property ("Check that strongClosure is closed (condition 2 of 3 for strong closure)") {
    forAll(GenSmallEvenInt) {
      (d: Int) =>
      forAll (GenFunMatrix(d)) {
        (m : FunMatrix[Double]) =>
        BagnaraStrongClosure.strongClosure(m) match {
          case None => false
          case Some(c) =>
            forAll (
              Gen.zip(
                Gen.choose(0,c.dimension.dim - 1),
                Gen.choose(0,c.dimension.dim - 1),
                Gen.choose(0,c.dimension.dim - 1))) {
              case (i, j, k) =>
                if (i == j)
                  c(i,j) == 0
                else
                  // No need to special case i = k or j = k
                  // it reduces to e.g. for i = k c(i,j) <= 0 + c(i,j)
                  c(i,j) <= c(i, k) + c(k, j)
            }
        }
      }
    }
  }

  property ("Check that for strongClosure m_ij <= (m_{i, bari}+m_{barj, j})/2 holds (condition 3 of 3 for strong closure)") {
    forAll(GenSmallEvenInt) {
      (d: Int) =>
      forAll (GenFunMatrix(d)) {
        (m : FunMatrix[Double]) =>
        BagnaraStrongClosure.strongClosure(m) match {
          case None => false
          case Some(c) =>
            forAll (
              Gen.zip(
                Gen.choose(0,c.dimension.dim - 1),
                Gen.choose(0,c.dimension.dim - 1))) {
              case (i, j) => {
                val ibar = signed(i)
                val jbar = signed(j)
                c(i,j) <= (c(i, ibar) + c(jbar, j)) / 2
              }
            }
        }
      }
    }
  }

  // Commented out as FunDBM is too inefficient for this to terminate in a reasonable time.
  // property ("Check that strongClosure produces legal values") {
  //   // i.e. this mainly means some NaNs we had to hunt down
  //   forAll(GenSmallEvenInt) {
  //     (d: Int) =>
  //     forAll (GenFunMatrix(d)) {
  //       case (m : FunMatrix[Double]) =>
  //         BagnaraStrongClosure.strongClosure(m) match {
  //           case None =>
  //             false
  //           case Some(c) =>
  //             checkIsLegal(c)
  //         }
  //     }
  //   }
  // }

  // Commented out as FunDBM is too inefficient for this to terminate in a reasonable time.
  // property ("Check that toInterval yields a valid interval") {
  //   forAll(GenSmallInt) {
  //     (d: Int) =>
  //     forAll(GenClosedFunDBMOrTop(d)) {
  //       case dbm : FunDBM[Closed, Double] =>
  //         {
  //           val o = new AbstractOctagon(dbm, oct, e)
  //           o.toInterval <= box.top(dbm.noOfVariables.count)
  //           o.toInterval >= box.bottom(dbm.noOfVariables.count)
  //         }
  //     }
  //   }
  // }

  // property ("Check that linearAssignment yields legal values") {
  //   forAll(GenSmallInt) {
  //     (d: Int) =>
  //     forAll(GenClosedFunDBMOrTop(d)) {
  //       case dbm : FunDBM[Closed, Double] =>
  //         {
  //           val o = new AbstractOctagon(dbm, oct, e)
  //           forAll(GenLf(o.dimension)) {
  //             case lf : LinearForm =>
  //               forAll(Gen.choose(0, o.dimension - 1)) {
  //                 case vi : Int =>
  //                   {
  //                     val ass = o.linearAssignment(vi, lf)
  //                     ass <= AbstractOctagon(e.topDBM[Double](VarCount(o.dimension)), oct, e)
  //                     ass.dbm match {
  //                       case dbm : ClosedFunDBM[Double] => checkIsLegal(dbm.m)
  //                       case b : BottomFunDBM[Double] => true
  //                       case _ => false
  //                     }
  //                   }
  //               }
  //           }
  //         }
  //     }
  //   }
  // }

  property ("Check that linearAssignment is sound, i.e. <= interval assignment") {
    forAll(GenSmallInt) {
      (d: Int) =>
      forAll(GenClosedFunDBMOrTop(d)) {
        case dbm : FunDBM[Closed, Double] =>
          {
            val o = new AbstractOctagon(dbm, oct, e)
            forAll(GenLf(o.dimension)) {
              case lf : LinearForm =>
                forAll(Gen.choose(0, o.dimension - 1)) {
                  case vi : Int =>
                    {
                      o.linearAssignment(vi, lf).toInterval <= o.toInterval.linearAssignment(vi, lf)
                    }
                }
            }
          }
      }
    }
  }

  property ("Check that forall X, Y : AbstractOctagon, (X widening Y) >= X, Y (condition 1 of 2 for soundness of widening)") {
    forAll(GenSmallInt) {
      (d: Int) =>
      forAll(GenClosedFunDBMOrTop(d)) {
        case dbmx : FunDBM[Closed, Double] =>
          {
            val x = new AbstractOctagon(dbmx, oct, e)
            forAll(GenClosedFunDBMOrTop(d)) {
              case dbmy : FunDBM[Closed, Double] =>
                {
                  val y = new AbstractOctagon(dbmy, oct, e)
                  (x widening y) >= x
                  (x widening y) >= y
                }
            }
          }
      }
    }
  }
}
