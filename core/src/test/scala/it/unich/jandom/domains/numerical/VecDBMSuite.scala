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
import org.scalatest.FunSuite
import variables._
import spire.math.Rational
import spire.math.RationalAlgebra
import it.unich.jandom.utils.numberext.RationalExt
import org.scalacheck.Arbitrary
import org.scalacheck.Arbitrary.arbitrary
import org.scalacheck.Gen

class VecDBMSuite extends FunSuite {
  val VecDBMInstance = (new DBMInstance[VecMatrix]()(VecMatrixMatrixInstance.vecMatrixIsMatrix))
  type VecDBM[B,C] = (DBM[VecMatrix, B, C])
  type ClosedVecDBM[A] = ClosedDBM[VecMatrix, A]
  type NonClosedVecDBM[A] = ClosedDBM[VecMatrix, A]
  val e = VecDBMInstance.funDBM
  val em =  VecMatrixMatrixInstance.vecMatrixIsMatrix
  val inf = RationalExt.PositiveInfinity
  val r = new RationalAlgebra()

  val box =  BoxRationalDomain()
  val Utils = new it.unich.jandom.domains.numerical.octagon.testutils.Utils(box)
  import Utils._
  import InfField.ifieldRationalExt
  val oct = new OctagonDomain[VecDBM, RationalExt, BoxRationalDomain](e, box)
  def id (x: RationalExt): RationalExt = x
  val dbmutil = new DBMUtils[RationalExt]()

  implicit def arbRational: Arbitrary[Rational] =
    Arbitrary {
      for {
        n <- arbitrary[Int]
        d <- arbitrary[Int]
      } yield(r.fromInt(n)) / r.fromInt(Math.max(1,Math.abs(d))) // max(1,d) is a hackish way to avoid n/0
    }


  implicit def arbBox : Arbitrary[box.Property] =
    Arbitrary {
      for {
        n <- Gen.choose(1,20)
        pairs : Array[(RationalExt, RationalExt)] <- Gen.containerOfN[Array, (RationalExt, RationalExt)](n, GenOrderedPair)
      } yield (new box.Property(pairs.unzip._1, pairs.unzip._2, false))
    }

  type DOM = OctagonDomain[FunDBM, RationalExt, BoxRationalDomain]



  test ("VecMatrix sanity check") {
    def f (i : Int, j : Int) : RationalExt = ((i + j) % 2)
    val m = em.make(f, Dimension(4))
    assert(m(0,0) == 0)
    assert(m(0,1) == 1)
    assert(m(3,3) == 0)
    def g (i : Int, j : Int) : RationalExt = ((i + j + 1) % 2)
    val n = m.update(g)
    assert(n(0,0) == 1)
    assert(n(0,1) == 0)
    def h(a : RationalExt, b : RationalExt) = a + b
    val o = m.combine(n, h)
    assert(o(0,0) == 1)
    assert(o(0,1) == 1)
    assert(o.toList == List(1,1,1,1,
      1,1,1,1,
      1,1,1,1,
      1,1,1,1))
  }

  test ("VecMatrixIsMatrix sanity check") {
    def f (i : Int, j : Int) : RationalExt = ((i + j) % 2)
    def g (i : Int, j : Int) : RationalExt = ((i + j + 1) % 2)
    def h(a : RationalExt, b : RationalExt) = a + b
    def i (x : RationalExt, z : => RationalExt) = (x + z)
    val m = em.make[RationalExt](f, Dimension(4))
    val n = em.make[RationalExt](g, Dimension(4))
    val fourfours = (em.pure[RationalExt](Dimension(2), 4))
    assert((em.combine(h)(m,n))(0,0) == 1)
    assert(fourfours.toList == List(4,4,4,4))
    assert(em.foldRight(fourfours, RationalExt(0))(i) == RationalExt(16))
  }

  test ("strongClosure sanity check") {
    // O4 and O4* from www.srl.inf.ethz.ch/pa2015/Octagon.pdf
    val o4 : Array[Array[RationalExt]]= Array(
      Array(0, -2,  0, -2,  inf,  inf),
      Array(2,  0,  2,  0,  inf,    0),
      Array(0, -2,  0, -2,  inf,  inf),
      Array(2,  0,  2,  0,  inf,  inf),
      Array(0, inf, inf, inf, 0, inf),
      Array(inf, inf, inf, inf, inf, 0)
    )

    val o4dbm : ClosedVecDBM[RationalExt] = ClosedDBM[VecMatrix, RationalExt](
      em.make[RationalExt]((i : Int, j : Int) => o4(i)(j), Dimension(6))
    )(em)

    val expectedo4star = Array[Array[RationalExt]](
      Array(0, -2, 0, -2, inf, -2),
      Array(2,  0, 2,  0, inf,  0),
      Array(0, -2, 0, -2, inf, -2),
      Array(2,  0, 2,  0, inf,  0),
      Array(0, -2, 0, -2,   0, -2),
      Array(inf, inf, inf, inf, inf, 0)
    )

    val o4star = e.strongClosure(o4dbm)

    o4star match {
      case dbm : ClosedVecDBM[RationalExt] => assert(dbm.m == em.make[RationalExt]((i,j) => expectedo4star(i)(j), Dimension(6)))
      case _ => assert(false)
    }
  }


  test ("T {v0 <- 2}.toInterval == [2,2]") {
    val a = AbstractOctagon(e.topDBM[RationalExt](VarCount(1)), oct, box, e)
    val c = LinearForm.c(2)
    val b = a.linearAssignment(0, c)
    assert(b.toInterval.high.head == 2)
    assert(b.toInterval.low.head == 2)
  }

  test ("Union of [1,1], [2, 2] == [1, 2]") {
    val a = AbstractOctagon(e.topDBM[RationalExt](VarCount(1)), oct, box, e)
    val c2 = LinearForm.c(1)
    val c1 = LinearForm.c(2)
    val b1 = a.linearAssignment(0, c1)
    val b2 = a.linearAssignment(0, c2)
    val union = b1 union b2
    assert(union.toInterval.isEmpty == false)
    assert(union.toInterval.high.head == 2)
    assert(union.toInterval.low.head == 1)
  }

  test ("Intersection of [1,1], [2,2] is empty") {
    val a = AbstractOctagon(e.topDBM[RationalExt](VarCount(1)), oct, box, e)
    val c2 = LinearForm.c(1)
    val c1 = LinearForm.c(2)
    val b1 = a.linearAssignment(0, c1)
    val b2 = a.linearAssignment(0, c2)
    val intersection = b1 intersection b2
    assert(intersection.toInterval.isEmpty == true)
  }

  test ("[1,2] <= [0,3]") {
    val a = AbstractOctagon(e.topDBM[RationalExt](VarCount(1)), oct, box, e)
    val c1 = LinearForm.c(1);  val c2 = LinearForm.c(2)
    val c3 = LinearForm.c(0);  val c4 = LinearForm.c(3)
    val b1 = a.linearAssignment(0, c1)
    val b2 = a.linearAssignment(0, c2)
    val i1 = b1 union b2 // i1 = [1,2]
    val b3 = a.linearAssignment(0, c3)
    val b4 = a.linearAssignment(0, c4)
    val i2 = b3 union b4 // i2 = [0,3]
    assert (i1 <= i2)
  }

  test ("Sanity check for decodeTest") {
    AbstractOctagon.decodeTest(
      new DenseLinearForm(Seq(0,0,0,0))
        // TODO Maybe this case is a special case of one of the others?
        // If so which one?
    ) match {
      case Fallback() => assert(true)
      case (other) => assert(false, "Got " + other)
    }

    AbstractOctagon.decodeTest(
      new DenseLinearForm(Seq(1,1,0,0))
    ) match {
      case Case1Test(vl, c) => assert(vl == 0 & c == c)
      case (other) => assert(false, "Got " + other)
    }

    AbstractOctagon.decodeTest(
      new DenseLinearForm(Seq(4,1,0,0))
    ) match {
      case Case1Test(vl, c) => assert(vl == 0 & c == c)
      case (other) => assert(false, "Got " + other)
    }

    AbstractOctagon.decodeTest(
      new DenseLinearForm(Seq(0,1,0,0))
    ) match {
      case Case1Test(vl, c) => assert(vl == 0 & c == c)
      case (other) => assert(false, "Got " + other)
    }

    AbstractOctagon.decodeTest(
      new DenseLinearForm(Seq(0,123,0,0))
    ) match {
      case Fallback() => assert(true)
      case (other) => assert(false, "Got " + other)
    }

    AbstractOctagon.decodeTest(
      new DenseLinearForm(Seq(1,0,-1,0))
    ) match {
      case Case2Test(vl, c) => assert(vl == 1 & c == 1)
      case (other) => assert(false, "Got " + other)
    }

    AbstractOctagon.decodeTest(
      new DenseLinearForm(Seq(1,1,0,-1))
    ) match {
      case Case3Test(vl, vk, c) => assert(vl == 0 & vk == 2 & c == 1)
      case (other) => assert(false, "Got " + other)
    }

    AbstractOctagon.decodeTest(
      new DenseLinearForm(Seq(1,1,0,1))
    ) match {
      case Case4Test(vl, vk, c) => assert(
        (vl == 0 & vk == 2 & c == 1)
          |
          (vl == 2 & vk == 0 & c == 1)
      )
      case (other) => assert(false, "Got " + other)
    }

    AbstractOctagon.decodeTest(
      new DenseLinearForm(Seq(1,-1,0,-1))
    ) match {
      case Case5Test(vl, vk, c) => assert(vl == 0 & vk == 2 & c == 1)
      case (other) => assert(false, "Got " + other)
    }

    AbstractOctagon.decodeTest(
      new DenseLinearForm(Seq(1,-1,1,-1))
    ) match {
      case Fallback() => assert(true)
      case (other) => assert(false, "Got " + other)
    }
  }

  test ("Sanity check for singleConstantExactAssignment") {
    val v = VarIndex(0)
    val c = 1
    val m = VecDBMInstance.funDBM.topDBM[RationalExt](VarCount(1))
    val r =
      e.incrementalClosure(v)(dbmutil.singleConstantExactAssignment(v, c)(m, e).elem)
    // See Vechev p8
    assert(r.innerMatrix.get(0,0) == 0)
    assert(r.innerMatrix.get(0,1) >= -2)
    assert(r.innerMatrix.get(1,1) == 0)
    assert(r.innerMatrix.get(1,0) >= 2)
  }



  test ("Sanity check for fallbackUpdate") {
    // CAVEAT: Super trivial example, not nearly enough to trust the thing
    val topoct = AbstractOctagon(e.topDBM[RationalExt](VarCount(3)), oct, box, e)
    val assigned =
      topoct.linearAssignment(0, LinearForm.c(1))
        .linearAssignment(1, LinearForm.c(2))
        .linearAssignment(2, LinearForm.c(3))
    // ie assign (1,2,3) to (v1, v2, v3)
    // now we assign Oct -= (v1 + v2 + v3), i.e we have
    // e = lf(0,1,1,1)
    val expectedassignment : Array[Array[RationalExt]] = Array(
      Array(0,   -2,   1, -3,   2, -4),
      Array(2,    0,   3, -1,   4, -2),
      Array(-1,  -3,   0, -4,   1, -5),
      Array( 3,   1,   4,  0,   5, -1),
      Array(-2,  -4,  -1, -5,   0, -6),
      Array( 4,   2,   5,  1,   6,  0)
    )

    assigned.dbm match {
      case dbm : ClosedVecDBM[RationalExt] =>  assert(dbm.m == em.make[RationalExt]((i,j) => expectedassignment(i)(j), Dimension(6)))
      case (other) => assert(false, "Got " + other)
    }

    val updated = assigned.fallbackUpdate(
      new DenseLinearForm(Seq(0,1,1,1))
    )

    val expected : Array[Array[RationalExt]] = Array(
      Array(0,  -10,  -5, -3,  -4, -4),
      Array(2,    0,  -9, -5, -10, -4),
      Array(-7,  -3,   0, -8,  -5, -5),
      Array(-9,  -7,   4,  0, -11, -5),
      Array(-8,  -4,  -7, -5,   0, -6),
      Array(-10, -8, -11, -7,   6,  0)
    )

    def peek[S, A](m: VecDBM[S, A]): Unit = m match {
      case dbm : ClosedDBM[VecMatrix, RationalExt] => assert(dbm.m == em.make[RationalExt]((i,j) => expected(i)(j), Dimension(6)))
      case dbm : NonClosedDBM[VecMatrix, RationalExt] => assert(dbm.m == em.make[RationalExt]((i,j) => expected(i)(j), Dimension(6)))
      case _ => assert(false, m)
    }
    peek(updated.elem)
  }

  test ("VecDBM-based AbstractOctagon matches example from www.srl.inf.ethz.ch/pa2015/Octagon.pdf p. 24") {
    /*
     * x <- 1
     * y <- x
     * while (x <= m)
     *   x <- x + 1
     *   y <- y + x
     * assert(y >= m)
     */
    val o1 = AbstractOctagon(e.topDBM[RationalExt](VarCount(3)), oct, box, e)

    // We choose (v0, v1, v2) = (x, y, m)
    // x <- 1
    val o2star = o1.linearAssignment(0, LinearForm.c(1))
    val expectedo2star = Array[Array[RationalExt]](
      Array(0, -2, inf, inf, inf, inf),
      Array(2,  0, inf,  inf, inf, inf),
      Array(inf, inf, 0, inf, inf, inf),
      Array(inf,  inf, inf,  0, inf, inf),
      Array(inf, inf, inf, inf, 0, inf),
      Array(inf, inf, inf, inf, inf, 0)
    )
    o2star.dbm match {
      case dbm : ClosedDBM[VecMatrix, RationalExt] => assert(dbm.m == em.make[RationalExt]((i: Int, j:Int) => expectedo2star(i)(j), Dimension(6)))
      case _ => assert(false)
    }

    // y <- x
    val o3star = o2star.linearAssignment(1, new DenseLinearForm(Seq(0,1,0,0)))
    val expectedo3star = Array[Array[RationalExt]](
      Array(0, -2, 0, -2, inf, inf),
      Array(2,  0, 2,  0, inf, inf),
      Array(0, -2, 0, -2, inf, inf),
      Array(2,  0, 2,  0, inf, inf),
      Array(inf, inf, inf, inf, 0, inf),
      Array(inf, inf, inf, inf, inf, 0)
    )

    o3star.dbm match {
      case dbm : ClosedDBM[VecMatrix, RationalExt] => assert(dbm.m == em.make[RationalExt]((i: Int, j:Int) => expectedo3star(i)(j), Dimension(6)))
      case _ => assert(false)
    }

    // (x <= m) ==> (x - m <= 0) == (v0 - v2 <= 0)

    val form = AbstractOctagon.decodeTest(new DenseLinearForm(Seq(0,1,0,-1)))
    form match {
      case Case3Test(vi, vj, c) => assert(vi == 0 & vj == 2)
      case (other) => assert(false, "Got " + other)
    }

    // First we check that the actual computation is correct...
    val o4 = o3star.linearInequalityEx(new DenseLinearForm(Seq(0,1,0,-1)))
    val expectedo4 = Array[Array[RationalExt]](
      Array(0, -2, 0, -2, inf, inf),
      Array(2,  0, 2,  0, inf, 0),
      Array(0, -2, 0, -2, inf, inf),
      Array(2,  0, 2,  0, inf, inf),
      Array(0, inf, inf, inf, 0, inf),
      Array(inf, inf, inf, inf, inf, 0)
    )

    def peek[S, A](m: VecDBM[S, A]): Unit = m match {
      case dbm : ClosedDBM[VecMatrix, RationalExt] => assert(dbm.m == em.make[RationalExt]((i,j) => expectedo4(i)(j), Dimension(6)))
      case dbm : NonClosedDBM[VecMatrix, RationalExt] => assert(dbm.m == em.make[RationalExt]((i,j) => expectedo4(i)(j), Dimension(6)))
      case _ => assert(false)
    }
    peek(o4.elem)

    // Then we check that linearEquality takes care of the closure
    val o4star = o3star.linearInequality(new DenseLinearForm(Seq(0,1,0,-1)))
    val expectedo4star = Array[Array[RationalExt]](
      Array(0, -2, 0, -2, inf, -2),
      Array(2,  0, 2,  0, inf, 0),
      Array(0, -2, 0, -2, inf, -2),
      Array(2,  0, 2,  0, inf, 0),
      Array(0, -2, 0, -2, 0, -2),
      Array(inf, inf, inf, inf, inf, 0)
    )
    o4star.dbm match {
      case dbm : ClosedDBM[VecMatrix, RationalExt] => assert(dbm.m == em.make[RationalExt]((i,j) => expectedo4star(i)(j), Dimension(6)))
      case _ => assert(false)
    }

    // If it's really closed it should be also a fixed point for the closure op?
    assert(e.strongClosure(o4star.dbm) == o4star.dbm)

    // x <- x + 1
    val o5star = o4star.linearAssignment(0, new DenseLinearForm(Seq(1,1,0,0)))
    val expectedo5star : Array[Array[RationalExt]] = Array(
      Array(0,  -4,  -1,  -3, inf, -3),
      Array(  4,   0,   3,   1, inf,  1),
      Array(  1,  -3,   0,  -2, inf, -2),
      Array(  3,  -1,   2,   0, inf,  0),
      Array(  1,  -3,   0,  -2,   0, -2),
      Array(inf, inf, inf, inf, inf,  0)
    )
    o5star.dbm match {
      case dbm : ClosedDBM[VecMatrix, RationalExt] => assert(dbm.m == em.make[RationalExt]((i,j) => expectedo5star(i)(j), Dimension(6)))
      case _ => assert(false)
    }

    // y <- y + x
    val o6star = o5star.linearAssignment(1, new DenseLinearForm(Seq(0,1,1,0)))
    val expectedo6star : Array[Array[RationalExt]] = Array(
      Array(  0,  -4,   1,  -5, inf, -5),
      Array(  4,   0,   5,  -1, inf, -1),
      Array( -1,  -5,   0,  -6, inf, -6),
      Array(  5,   1,   6,   0, inf,  0),
      Array(  0,  -4,   1,  -5,   0, -4),
      Array(inf, inf, inf, inf, inf,  0)
    )

    // o6star.dbm match {
    //   case dbm : ClosedDBM[VecMatrix, RationalExt] => assert(dbm.m == em.make[RationalExt]((i,j) => expectedo6star(i)(j), 6))
    //   case _ => assert(false)
    // }
    //
    // TODO This doesn't hold, probably because Vechev uses a different approximation in the slides? (this is one of the "interval fallback" cases)

    val expectedo3star2 = Array[Array[RationalExt]](
      Array(0, -2, 1, -2, inf, inf),
      Array(4,  0, 5,  0, inf, inf),
      Array(0, -2, 0, -2, inf, inf),
      Array(5,  1, 6,  0, inf, inf),
      Array(inf, inf, inf, inf, 0, inf),
      Array(inf, inf, inf, inf, inf, 0)
    )

    o6star.union(o3star).dbm match {
      case dbm : ClosedDBM[VecMatrix, RationalExt] => assert(dbm.m == em.make[RationalExt]((i,j) => expectedo3star2(i)(j), Dimension(6)))
      case _ => assert(false)
    }
  }
}
