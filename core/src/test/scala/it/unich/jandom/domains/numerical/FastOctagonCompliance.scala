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

class FastOctagonCompliance extends PropSpec with PropertyChecks {
  val mev = MEvidence(
    HalfMatrixDenseSparseInstance.halfMatrixDecomposableInstance,
    HalfMatrixDenseSparseInstance.halfMatrixDenseSparseInstance,
    HalfMatrixDenseSparseInstance.halfMatrixDecomposableInstance)
  val fe = CFDBMInstance.instance(mev)
  type FastDBM[S, A] = CFastDBM[HalfMatrix, HalfSubMatrix, S, A]
  val box =  BoxRationalDomain()
  val foct = new OctagonDomain[FastDBM, RationalExt, BoxRationalDomain](fe, box)
  val Utils = new it.unich.jandom.domains.numerical.octagon.testutils.Utils(box)
  def id (x: RationalExt): RationalExt = x
  import Utils._

  implicit def arbBox : Arbitrary[box.Property] =
    Arbitrary {
      for {
        n <- Gen.choose(1,20)
        pairs : Array[(RationalExt, RationalExt)] <- Gen.containerOfN[Array, (RationalExt, RationalExt)](n, GenOrderedPair)
      } yield (new box.Property(pairs.unzip._1, pairs.unzip._2, false))
    }

  val VecDBMInstance = (new DBMInstance[VecMatrix]()(VecMatrixMatrixInstance.vecMatrixIsMatrix))
  val e = VecDBMInstance.funDBM
  val oct = new OctagonDomain[VecDBM, RationalExt, BoxRationalDomain](e, box)
  type VecDBM[B,C] = (DBM[VecMatrix, B, C])

  def GenVarAndLf(n: Int) : Gen[(LinearForm, Int)] = {
    require (n > 0)
    Gen.zip (GenLf(n), Gen.choose(0, n-1))
  }

  def GenVarAndLfAndIneq(n: Int) : Gen[(LinearForm, Int, IsInequality)] = {
    require (n > 0)
    Gen.zip (GenLf(n), Gen.choose(0, n-1), arbitrary[Boolean])
  }

  def GenListOfVarAndLf(n: Int) : Gen[Seq[(LinearForm, Int)]] = for {
    forms <- GenSmallInt
    list <- Gen.listOfN(forms, GenVarAndLf(n))
  } yield list


  type IsInequality = Boolean
  def GenListOfVarAndLfAndIneq(n: Int) : Gen[Seq[(LinearForm, Int, IsInequality)]] = for {
    forms <- GenSmallInt
    list <- Gen.listOfN(forms, GenVarAndLfAndIneq(n))
  } yield list


  type FDOM = OctagonDomain[FastDBM, RationalExt, BoxRationalDomain]
  type DOM = OctagonDomain[VecDBM, RationalExt, BoxRationalDomain]

  //////////////////////////////////////////////////////////////////////////////
  // Begin properties
  //////////////////////////////////////////////////////////////////////////////

  //  def noShrink[T](gen: Gen[T]): Gen[NoShrinkWrapper[T]] = gen.map(NoShrinkWrapper.apply)
  //  case class NoShrinkWrapper[T](value: T)
  import org.scalacheck.Shrink
  implicit val noShrink: Shrink[Int] = Shrink.shrinkAny

  def compare(
    octf:  AbstractOctagon[FDOM, FastDBM,  RationalExt, BoxRationalDomain],
    octmine: AbstractOctagon[DOM, VecDBM,  RationalExt, BoxRationalDomain]) : Boolean = {
    if (octf.dimension != octmine.dimension) false
    else
      (0 until octf.dimension).forall{
        (i: Int) =>
        (0 until octf.dimension).forall{
          (j: Int) => fe.get(i,j)(octf.dbm) == e.get(i,j)(octmine.dbm)
        }
      }
  }

  type FastOct = AbstractOctagon[FDOM, FastDBM,  RationalExt, BoxRationalDomain]
  type MineOct = AbstractOctagon[DOM, VecDBM,  RationalExt, BoxRationalDomain]

  property ("Fast octagons match Mine octagons in their response to linear assignment") {
      forAll(GenSmallInt) {
        (d: Int) => {
          forAll(GenListOfVarAndLf(d)) {
            (s: Seq[(LinearForm, Int)]) => {
              val topf =  AbstractOctagon[FDOM, FastDBM,  RationalExt, BoxRationalDomain](fe.topDBM[RationalExt](VarCount(d)), foct, box, fe)
              val topmine = AbstractOctagon[DOM, VecDBM,  RationalExt, BoxRationalDomain](VecDBMInstance.funDBM.topDBM[RationalExt](VarCount(d)), oct, box, e)

              def f(s: Seq[(LinearForm, Int)],
                octf:  FastOct,
                octmine: MineOct) : Boolean = {
                if (s.size == 0)
                  true
                else
                  if (compare(octf, octmine)) {
                    f(
                      s.tail,
                      octf.linearAssignment(s.head._2, s.head._1),
                      octmine.linearAssignment(s.head._2, s.head._1)
                    )
                  }
                  else
                    false
              }
              f(s, topf, topmine)
            }
          }
        }
      }
  }



  property ("Fast octagons match Mine octagons in their response to linear assignment AND inequality") {
      forAll(GenSmallInt) {
        (d: Int) => {
          forAll(GenListOfVarAndLfAndIneq(d)) {
            (s: Seq[(LinearForm, Int, IsInequality)]) => {
              val topf =  AbstractOctagon[FDOM, FastDBM,  RationalExt, BoxRationalDomain](fe.topDBM[RationalExt](VarCount(d)), foct, box, fe)
              val topmine = AbstractOctagon[DOM, VecDBM,  RationalExt, BoxRationalDomain](VecDBMInstance.funDBM.topDBM[RationalExt](VarCount(d)), oct, box, e)

              def f(s: Seq[(LinearForm, Int, IsInequality)],
                octf:  FastOct,
                octmine: MineOct) : Boolean = {
                if (s.size == 0)
                  true
                else
                  if (compare(octf, octmine)) {
                    if (s.head._3) {
                      f(
                        s.tail,
                        octf.linearInequality(s.head._1),
                        octmine.linearInequality(s.head._1)
                      )
                    } else {
                      f(
                        s.tail,
                        octf.linearAssignment(s.head._2, s.head._1),
                        octmine.linearAssignment(s.head._2, s.head._1)
                      )
                    }
                  }
                  else
                    false
              }


              f(s, topf, topmine)
            }
          }
        }
      }
  }

  property ("Fast octagons match Mine octagons in their response to union") {
      forAll(GenSmallInt) {
        (d: Int) => {
          forAll(GenListOfVarAndLfAndIneq(d)) {
            (s: Seq[(LinearForm, Int, IsInequality)]) => {
              val topf =  AbstractOctagon[FDOM, FastDBM,  RationalExt, BoxRationalDomain](fe.topDBM[RationalExt](VarCount(d)), foct, box, fe)
              val topmine = AbstractOctagon[DOM, VecDBM,  RationalExt, BoxRationalDomain](VecDBMInstance.funDBM.topDBM[RationalExt](VarCount(d)), oct, box, e)

              def f(s: Seq[(LinearForm, Int, IsInequality)],
                octf:  FastOct,
                octmine: MineOct) : (FastOct, MineOct) = {
                if (s.size == 0)
                  (octf, octmine)
                else {
                  assert (compare(octf, octmine))
                  if (s.head._3) {
                    f(
                      s.tail,
                      octf.linearInequality(s.head._1),
                      octmine.linearInequality(s.head._1)
                    )
                  } else {
                    f(
                      s.tail,
                      octf.linearAssignment(s.head._2, s.head._1),
                      octmine.linearAssignment(s.head._2, s.head._1)
                    )
                  }
                }
              }

              val (octf1, octmine1) = f(s, topf, topmine)
              val (octf2, octmine2) = f(s, topf, topmine)
              assert(compare(octf1.union(octf2), octmine1.union(octmine2)))
            }
          }
        }
      }
  }

  property ("Fast octagons match Mine octagons in their response to intersection") {
      forAll(GenSmallInt) {
        (d: Int) => {
          forAll(GenListOfVarAndLfAndIneq(d)) {
            (s: Seq[(LinearForm, Int, IsInequality)]) => {
              val topf =  AbstractOctagon[FDOM, FastDBM,  RationalExt, BoxRationalDomain](fe.topDBM[RationalExt](VarCount(d)), foct, box, fe)
              val topmine = AbstractOctagon[DOM, VecDBM,  RationalExt, BoxRationalDomain](VecDBMInstance.funDBM.topDBM[RationalExt](VarCount(d)), oct, box, e)

              def f(s: Seq[(LinearForm, Int, IsInequality)],
                octf:  FastOct,
                octmine: MineOct) : (FastOct, MineOct) = {
                if (s.size == 0)
                  (octf, octmine)
                else {
                  assert (compare(octf, octmine))
                  if (s.head._3) {
                    f(
                      s.tail,
                      octf.linearInequality(s.head._1),
                      octmine.linearInequality(s.head._1)
                    )
                  } else {
                    f(
                      s.tail,
                      octf.linearAssignment(s.head._2, s.head._1),
                      octmine.linearAssignment(s.head._2, s.head._1)
                    )
                  }
                }
              }

              val (octf1, octmine1) = f(s, topf, topmine)
              val (octf2, octmine2) = f(s, topf, topmine)
              assert(compare(octf1.intersection(octf2), octmine1.intersection(octmine2)))
            }
          }
        }
      }
  }
}
