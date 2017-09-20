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
import it.unich.jandom.domains.numerical._
import it.unich.jandom.domains.numerical.octagon._
import it.unich.jandom.domains.numerical.octagon.variables.VarCount
import InfField._
import org.scalatest.FunSuite
import spire.math.Rational
import it.unich.jandom.utils.numberext.RationalExt
import variables.{VarIndex, Dimension}

class FastOctagonSuite extends FunSuite {

  val mev = MEvidence(
    HalfMatrixDenseSparseInstance.halfMatrixDecomposableInstance,
    HalfMatrixDenseSparseInstance.halfMatrixDenseSparseInstance,
    HalfMatrixDenseSparseInstance.halfMatrixDecomposableInstance)

  val inf = Double.PositiveInfinity
  val e = CFDBMInstance.instance(mev)
  type FastDBM[S, A] = CFastDBM[HalfMatrix, HalfSubMatrix, S, A]
  val box =  BoxRationalDomain()
  val oct = new OctagonDomain[FastDBM, RationalExt, BoxRationalDomain](e, box)
  type DOM = OctagonDomain[FastDBM, RationalExt, BoxRationalDomain]

  test("Sanity check for pure") {
    assert(mev.dec.pure(VarCount(1), RationalExt.PositiveInfinity).vec.size == 4)
    assert(mev.dec.pure(VarCount(1), RationalExt.PositiveInfinity).vec.forall(_ == RationalExt.PositiveInfinity))
    assert(mev.dec.pure(VarCount(2), RationalExt.PositiveInfinity).vec.size == 12)
    assert(mev.dec.pure(VarCount(2), 0).vec.forall(_ == 0))
  }

  test("T {v0 <- 123}.toInterval == [123,123]") {
    val top = e.topDBM[RationalExt](VarCount(1))
    val t = AbstractOctagon[DOM, FastDBM, RationalExt, BoxRationalDomain](top, oct, box, e)
    assert(t.dbm == top)
    val a = t.linearAssignment(0, LinearForm.c(123))
    assert(a.toInterval.high.head == 123 & a.toInterval.low.head == 123)
  }

  // TODO Once found the cause, rename
  test("bug #697659 - failing requirement") {
    val top = e.topDBM[RationalExt](VarCount(2))
    val t = AbstractOctagon[DOM, FastDBM, RationalExt, BoxRationalDomain](top, oct, box, e)
    assert(t.dbm == top)
    val working = t.linearAssignment(0, DenseLinearForm(Seq(Rational(1)/Rational(1),Rational(1)/Rational(1))))
    val failing = t.linearAssignment(0, DenseLinearForm(Seq(Rational(1)/Rational(1),Rational(2)/Rational(1))))
  }

  // TODO Once found the cause, rename
  test("bug #697684 - failing requirement") {
    val top = e.topDBM[RationalExt](VarCount(2))
    val t = AbstractOctagon[DOM, FastDBM,  RationalExt, BoxRationalDomain](top, oct, box, e)
    assert(t.dbm == top)
    val foo = DenseLinearForm(Seq(Rational(0), Rational(1)))
    val bar = DenseLinearForm(Seq(Rational(-1)/Rational(1), Rational(-1)/Rational(1)))
    val first = t.linearAssignment(1, foo)
    val second = first.linearInequality(bar)
  }

  test("Sanity check for calculateComponents") {
    val top = e.topDBM[RationalExt](VarCount(4))
    val t = AbstractOctagon[DOM, FastDBM,  RationalExt, BoxRationalDomain](top, oct, box, e)
    assert(t.dbm == top)
    // At first, T_4 [v0 <- c] has 1 independent component
    val first = t.linearAssignment(0, DenseLinearForm(Seq(Rational(0))))
    first.dbm match {
      case CFast(m) => {
        val comp = FastDbmUtils.calculateComponents(m)(mev, ifieldRationalExt)
        assert(comp == List(List(VarIndex(0))))
      }
      case _ => assert(false)
    }
    // T_4 [v0 <- c][v3 <- c] has 1 independent component
    val second = first.linearAssignment(3, DenseLinearForm(Seq(Rational(1))))
    second.dbm match {
      case CFast(m) => {
        val comp = FastDbmUtils.calculateComponents(m)(mev, ifieldRationalExt)
        assert(comp.map(_.toSet).toSet == Set(Set(VarIndex(0), VarIndex(3))))
      }
      case _ => assert(false)
    }
    // T_4 [v0 <- c][v3 <- c][v3 <- v0] has 1 independent component again
    val third = second.linearAssignment(3, DenseLinearForm(Seq(Rational(0), Rational(1))))
    third.dbm match {
      case CFast(m) => {
        val comp = FastDbmUtils.calculateComponents(m)(mev, ifieldRationalExt)
        assert(comp.map(_.toSet).toSet == Set(Set(VarIndex(0), VarIndex(3))))
      }
      case _ => assert(false)
    }
  }

  test ("Check calculateComponents works with fig 3 from Singh, Pueschel, Vechev 2015") {
    val oo = ifieldRationalExt.infinity
    val a = Array[Array[RationalExt]](
      Array(0, oo),
      Array(oo, 0),
      Array(oo, oo, 0, 1),
      Array(oo, oo, 0, 0),
      Array(2,   1, oo, oo, 0, oo),
      Array(oo, oo, oo, oo, 2,  0),
      Array(oo, oo, oo, oo, oo, oo, 0, oo),
      Array(oo, oo, oo, oo, oo, oo, oo, 0),
      Array(oo, oo, oo, oo,  1, oo, oo, oo, 0,  4),
      Array(oo, oo, oo, oo,  2, oo, oo, oo, oo, 0)
    )

    val example = e.fromFun(Dimension(10), (i,j) => a(i)(j))

    example match {
      case CFast(m) => {
        val comp = FastDbmUtils.calculateComponents[HalfMatrix, HalfSubMatrix, RationalExt](m)(mev, ifieldRationalExt)
        assert(comp.map(_.toSet).toSet == Set(Set(VarIndex(0), VarIndex(2), VarIndex(4)), Set(VarIndex(1))))
      }
      case _ => assert(false)
    }
  }

  test ("Union of T, T == T") {
    val top = AbstractOctagon[DOM, FastDBM,  RationalExt, BoxRationalDomain](e.topDBM[RationalExt](VarCount(1)), oct, box, e)
    val union = top union top
    assert(union.toInterval.isEmpty == false)
    assert(union == top)
  }

  test ("Union of T, _|_ == _|_") {
    val top = AbstractOctagon[DOM, FastDBM,  RationalExt, BoxRationalDomain](e.topDBM[RationalExt](VarCount(1)), oct, box, e)
    val bot = AbstractOctagon[DOM, FastDBM,  RationalExt, BoxRationalDomain](e.bottomDBM[RationalExt](VarCount(1)), oct, box, e)
    val union = top union bot
    assert(union.toInterval.isEmpty == true)
    assert(union == bot)
  }

  test ("Union of [1,1], [2, 2] == [1, 2]") {
    val top = AbstractOctagon[DOM, FastDBM,  RationalExt, BoxRationalDomain](e.topDBM[RationalExt](VarCount(1)), oct, box, e)
    val c2 = LinearForm.c(1)
    val c1 = LinearForm.c(2)
    val b1 = top.linearAssignment(0, c1)
    val b2 = top.linearAssignment(0, c2)
    val union = b1 union b2
    assert(union.toInterval.isEmpty == false)
    assert(union.toInterval.high.head == 2)
    assert(union.toInterval.low.head == 1)
  }

  test ("Intersection of [1,1], [2,2] is empty") {
    val top = AbstractOctagon[DOM, FastDBM,  RationalExt, BoxRationalDomain](e.topDBM[RationalExt](VarCount(1)), oct, box, e)
    val c2 = LinearForm.c(1)
    val c1 = LinearForm.c(2)
    val b1 = top.linearAssignment(0, c1)
    val b2 = top.linearAssignment(0, c2)
    val intersection = b1 intersection b2
    assert(intersection.toInterval.isEmpty == true)
  }

  test ("Diagonal of union is 0") {
    val top = AbstractOctagon[DOM, FastDBM,  RationalExt, BoxRationalDomain](e.topDBM[RationalExt](VarCount(4)), oct, box, e)
    val b1 = top.linearAssignment(0, 123)
      .linearAssignment(1, 123)
      .linearAssignment(2, 123)
      .linearAssignment(3, 123)

    val b2 = top.linearAssignment(0, 456)
      .linearAssignment(1, 456)
      .linearAssignment(2, 456)
      .linearAssignment(3, 456)

    assert(e.get(0,0)(b1.dbm).get == 0)
    assert(e.get(1,1)(b1.dbm).get == 0)
    assert(e.get(2,2)(b1.dbm).get == 0)
    assert(e.get(3,3)(b1.dbm).get == 0)

    assert(e.get(0,0)(b2.dbm).get == 0)
    assert(e.get(1,1)(b2.dbm).get == 0)
    assert(e.get(2,2)(b2.dbm).get == 0)
    assert(e.get(3,3)(b2.dbm).get == 0)

    val union = b1 union b2
    assert(e.get(0,0)(union.dbm).get == 0)
    assert(e.get(1,1)(union.dbm).get == 0)
    assert(e.get(2,2)(union.dbm).get == 0)
    assert(e.get(3,3)(union.dbm).get == 0)
  }

  test ("T[v1 <- 2v0 - 1][v3 <- 2v2 - 1] U  T[v1 <- v0][v3 <- v2] has independen components [{v1, v0}, {v2, v3}}") {
    val top = e.topDBM[RationalExt](VarCount(4))
    val expectedComponents = Set(Set(VarIndex(0), VarIndex(1)), Set(VarIndex(2), VarIndex(3)))
    val t = AbstractOctagon[DOM, FastDBM,  RationalExt, BoxRationalDomain](top, oct, box, e)
    // T[v1 <- v0][v3 <- v2]
    val ass = t.linearAssignment(1, DenseLinearForm(Seq(1,1))).linearAssignment(3, DenseLinearForm(Seq(1,0,0,1)))
    ass.dbm match {
      case CFast(m) => {
        val comp = FastDbmUtils.calculateComponents(m)(mev, ifieldRationalExt)
        assert(comp.map(_.toSet).toSet == expectedComponents)
      }
      case _ => assert(false)
    }
    // T[v1 <- 2v0 - 1][v3 <- 2v2 - 1]
    val ass2 = t.linearAssignment(1, DenseLinearForm(Seq(-1,1))).linearAssignment(3, DenseLinearForm(Seq(-1,0,0,1)))
    ass2.dbm match {
      case CFast(m) => {
        val comp = FastDbmUtils.calculateComponents(m)(mev, ifieldRationalExt)
        assert(comp.map(_.toSet).toSet == expectedComponents)
      }
      case _ => assert(false)
    }
    val union = ass.union(ass2)
    union.dbm match {
      case CFast(m) => {
        val comp = FastDbmUtils.calculateComponents(m)(mev, ifieldRationalExt)
        assert(comp.map(_.toSet).toSet == expectedComponents)
      }
      case _ => assert(false)
    }
  }
}
