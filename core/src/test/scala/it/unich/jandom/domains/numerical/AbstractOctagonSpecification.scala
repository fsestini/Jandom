package it.unich.jandom.domains.numerical.octagon

import it.unich.jandom.domains.numerical.octagon._
import it.unich.jandom.domains.numerical._
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

trait AbstractOctagonSpecification extends PropSpec with PropertyChecks {

  type A
  type M[_,_]

  val utils = new DBMUtils[A]

  def GenFinite: Gen[A]
  def GenMatrix(dimension: Int, probInf: Int = 20): Gen[M[NonClosed, A]]

  def e: DifferenceBoundMatrix[M] { type PosetConstraint[X] = InfField[X] }

  type B <: BoxGenericDomain[A]
  val box: B
  def oct = new OctagonDomain[M, A, B](e, box)

  implicit def ifield: InfField[A]
  implicit def ratToN: (RationalExt => A)
  def toRat: (A => Rational)
  implicit def ct: ClassTag[A]
  // implicit def arbA: Arbitrary[A]
  // implicit def arbBox: Arbitrary[box.Property]

  def createInterval(low: Array[A], high: Array[A]): box.Property

  lazy val Utils = new AbstractUtils[A, B](box, ifield, GenFinite)
  import Utils._

  implicit def arbBox : Arbitrary[box.Property] =
    Arbitrary {
      for {
        n <- Gen.choose(1,20)
        pairs: Array[(A, A)]
          <- Gen.containerOfN[Array, (A, A)](n, GenOrderedPair)
      } yield (createInterval(pairs.unzip._1, pairs.unzip._2))
    }

  type DOM = OctagonDomain[M, A, B]

  def GenClosedMatrix(d: Int): Gen[M[Closed, A]] = for {
    m <- GenMatrix(d)
  } yield {
    // so sad that you have to spell out *everything*
    val c: M[Closed, A] = e.strongClosure(m)
    c
  }

  def GenClosedNonBottomMatrix(d: Int): Gen[M[Closed, A]] =
    GenClosedMatrix(d) suchThat (!e.isBottomDBM(_))

  def GenOctagon(d: Int): Gen[AbstractOctagon[DOM,M,A,B]] = for {
    cm <- GenClosedMatrix(d)
  } yield {
    // so sad that you have to spell out *everything*
    val o: AbstractOctagon[DOM,M,A,B] =
      new AbstractOctagon[DOM, M, A, B](cm, oct, box, e)
    o
  }



  def lt(x: A, y: A): Boolean = ifield.compare(x, y) == LT
  def gt(x: A, y: A): Boolean = ifield.compare(x, y) == GT
  def leq(x: A, y: A): Boolean =
    ifield.compare(x, y) == LT || ifield.compare(x, y) == EQ
  def eq(x: A, y: A): Boolean = ifield.compare(x, y) == EQ

  def coherent(d: Int, cc: M[Closed, A]) =
    forAll (Gen.zip(Gen.choose(0,d - 1), Gen.choose(0,d - 1))) {
      case (i, j) =>
        val ibar = signed(i)
        val jbar = signed(j)
          (e.get(jbar, ibar)(cc), e.get(i, j)(cc)) match {
          case (Some(x), Some(y)) => x == y
          case (None, None) => true
          case _ => false
        }
    }

  def closed(d: Int, closed: M[Closed, A]) =
    forAll (GenMatrix(d)) { (m : M[NonClosed, A]) =>
      val closed = e.strongClosure(m)
      forAll (Gen.zip(
        Gen.choose(0,d - 1), Gen.choose(0,d - 1), Gen.choose(0,d - 1))) {
        case (i, j, k) =>
          (e.get(i,j)(closed), e.get(i,k)(closed), e.get(k,j)(m)) match {
            case (Some(ij), Some(ik), Some(kj)) =>
              if (i == j) eq(ij, ifield.zero) else leq(ij, ifield.+(ik, kj))
            case (None, None, None) => true
            case _ => false
          }
      }
    }

  def condition3(d: Int, closed: M[Closed, A]) =
    forAll (Gen.zip(Gen.choose(0,d - 1), Gen.choose(0,d - 1))) {
      case (i, j) => {
        val ibar = signed(i)
        val jbar = signed(j)
          (e.get(i,j)(closed), e.get(i,ibar)(closed), e.get(jbar,j)(closed)) match {
          case (Some(ij), Some(iibar), Some(jbarj)) =>
            leq(ij, ifield.half(ifield.+(iibar, jbarj)))
          case (None, None, None) => true
          case _ => false
        }
      }
    }

  property("Can create from interval") {
    forAll { (b1: box.Property) =>
      val a = AbstractOctagon.fromInterval[DOM, M, A, B](b1, oct, box, e)
      true
    }
  }

  property ("T.{x <- C} toInterval is [C,C]") {
    forAll { (c: Int) =>
      val a = AbstractOctagon(e.topDBM[A](VarCount(1)), oct, box, e)
      val b = a.linearAssignment(0, LinearForm.c(c))
      (b.toInterval.high.head == c && b.toInterval.low.head == c)
    }
  }

  property("toInterval.fromInterval == id") {
    forAll { (b1: box.Property) =>
      AbstractOctagon
        .fromInterval[DOM, M, A, B](b1, oct, box, e)
        .toInterval == b1
    }
  }

  property ("T >= _|_") {
    forAll (GenSmallInt) { (d: Int) =>
      val topoct = AbstractOctagon(e.topDBM[A](VarCount(d)), oct, box, e)
      val botoct = AbstractOctagon(e.bottomDBM[A](VarCount(d)), oct, box, e)
      (topoct >= botoct)
    }
  }

  property ("T <> _|_ for different dim") {
    forAll (GenSmallInt) { (d: Int) =>
      forAll (GenSmallInt) { (delta: Int) =>
        val d1 = d + delta
        val topoct =  AbstractOctagon(e.topDBM[A](VarCount(d)), oct, box, e)
        val botoct = AbstractOctagon(e.bottomDBM[A](VarCount(d1)), oct, box, e)
        intercept[IllegalArgumentException]{topoct >= botoct}
        intercept[IllegalArgumentException]{botoct >= topoct}
      }
    }
  }

  property ("T <> T for different dim") {
    forAll (GenSmallInt) { (d: Int) =>
      forAll (GenSmallInt) { (delta: Int) =>
        val d1 = d + delta
        val topoct =  AbstractOctagon(e.topDBM[A](VarCount(d)), oct, box, e)
        val topoct2 = AbstractOctagon(e.topDBM[A](VarCount(d1)), oct, box, e)
        intercept[IllegalArgumentException]{topoct >= topoct2}
        intercept[IllegalArgumentException]{topoct2 >= topoct}
      }
    }
  }


  property ("T <- c results in [c,c]") {
    forAll { (c: Int) =>
      val a = AbstractOctagon(e.topDBM[A](VarCount(1)), oct, box, e)
      val b = a.linearAssignment(0, LinearForm.c(c))
      (b.toInterval.high.head == c && b.toInterval.low.head == c)
    }
  }

  property ("_|_ <> _|_ for different dim") {
    forAll (GenSmallInt) { (d: Int) =>
      forAll (GenSmallInt) { (delta: Int) =>
        val d1 = d + delta
        val botoct =  AbstractOctagon(e.bottomDBM[A](VarCount(d)), oct, box, e)
        val botoct2 = AbstractOctagon(e.bottomDBM[A](VarCount(d1)), oct, box, e)
        intercept[IllegalArgumentException]{botoct >= botoct2}
        intercept[IllegalArgumentException]{botoct2 >= botoct}
      }
    }
  }

  property ("T >= T") {
    forAll (GenSmallInt) { (d: Int) =>
      val topoct =  AbstractOctagon(e.topDBM[A](VarCount(d)), oct, box, e)
      val anotherTopoct = AbstractOctagon(e.bottomDBM[A](VarCount(d)), oct, box, e)
      (topoct >= anotherTopoct)
    }
  }

  property ("T U _|_ = T != _|_") {
    forAll (GenSmallInt) { (d: Int) =>
      val topoct = AbstractOctagon(e.topDBM[A](VarCount(d)), oct, box, e)
      val botoct = AbstractOctagon(e.bottomDBM[A](VarCount(d)), oct, box, e)
      ((topoct union botoct) == topoct && (topoct union botoct) != botoct)
    }
  }

  property ("T and _|_ intersection T = _|_") {
    forAll (GenSmallInt) { (d: Int) =>
      val topoct = AbstractOctagon(e.topDBM[A](VarCount(d)), oct, box, e)
      val botoct = AbstractOctagon(e.bottomDBM[A](VarCount(d)), oct, box, e)
      ((topoct intersection botoct) != topoct &&
        (topoct intersection botoct) == botoct)
    }
  }

  property ("Union of [C1,C1], [C2, C2] == [C1, C2] w/C1 < C2") {
    forAll(GenOrderedDistinctFinitePair) { (c : (A, A)) =>
      val a = AbstractOctagon(e.topDBM[A](VarCount(1)), oct, box, e)
      val b1 = a.linearAssignment(0, LinearForm.c(toRat(c._1)))
      val b2 = a.linearAssignment(0, LinearForm.c(toRat(c._2)))
      val union = b1 union b2
      assert(ifield.compare(c._1,c._2) == LT)
      (union.toInterval.isEmpty == false &&
        union.toInterval.high.head == c._2 &&
        union.toInterval.low.head == c._1)
    }
  }

  property ("Intersection of [C1,C1], [C2, C2] == _|_ w/C1 < C2") {
    forAll(GenOrderedDistinctFinitePair) { (c : (A, A)) =>
      val a = AbstractOctagon(e.topDBM[A](VarCount(1)), oct, box, e)
      val b1 = a.linearAssignment(0, LinearForm.c(toRat(c._1)))
      val b2 = a.linearAssignment(0, LinearForm.c(toRat(c._2)))
      val intersection = b1 intersection b2
      assert(ifield.compare(c._1, c._2) == LT)
      (intersection.isBottom == true &&
        intersection.toInterval.isEmpty == true)
    }
  }

  property ("[C1,C2] <= [C3<C1, C4>C2]") {
    forAll(GenOrderedDistinctFinitePair) { (c: (A, A)) =>
      forAll(GenOrderedDistinctFinitePair) { (d: (A, A)) =>
        val c1 = c._1
        val c2 = c._2
        assert(lt(c._1, c._2))
        assert(lt(d._1, d._2))
        val posC = (ifield.-(d._2, d._1)) // is certainly positive
        assert(gt(posC, ifield.zero))
        val c3 = ifield.-(c1, posC)
        val c4 = ifield.+(c2, posC)
        assert(lt(c3, c1) && lt(c1, c4) && gt(c4, c2))
        val a = AbstractOctagon(e.topDBM[A](VarCount(1)), oct, box, e)
        val b1 = a.linearAssignment(0, LinearForm.c(toRat(c1)))
        val b2 = a.linearAssignment(0, LinearForm.c(toRat(c2)))
        val b3 = a.linearAssignment(0, LinearForm.c(toRat(c3)))
        val b4 = a.linearAssignment(0, LinearForm.c(toRat(c4)))
        val i1 = b1 union b2 // i1 = [c1,c2]
        val i2 = b3 union b4 // i2 = [c3,c4]
        (i1 <= i2)
      }
    }
  }

  property ("Check that strongClosure is coherent (condition 1 of 3 for strong closure)") {
    forAll(GenSmallEvenInt) { (d: Int) =>
      forAll (GenClosedMatrix(d)) { (m : M[Closed, A]) =>
        coherent(d, m)
      }
    }
  }

  property ("Check that strongClosure is closed (condition 2 of 3 for strong closure)") {
    forAll(GenSmallEvenInt) { (d: Int) =>
      forAll (GenClosedMatrix(d)) { (m : M[Closed, A]) =>
        closed(d, m)
      }
    }
  }

  property ("Check that for strongClosure m_ij <= (m_{i, bari}+m_{barj, j})/2 holds (condition 3 of 3 for strong closure)") {
    forAll(GenSmallEvenInt) { (d: Int) =>
      forAll (GenClosedMatrix(d)) { (m : M[Closed, A]) =>
        condition3(d, m)
      }
    }
  }

  property ("Check that linearAssignment is sound, i.e. <= interval assignment") {
    forAll(GenSmallEvenInt) { (dd: Int) =>
      val d = dd + 2
      forAll(GenMatrix(d)) { (m : M[NonClosed, A]) =>
        val dbm = e.strongClosure(m)
        val o = new AbstractOctagon[DOM, M, A, B](dbm, oct, box, e)
        forAll(GenLf(o.dimension)) { (lf : LinearForm) =>
          forAll(Gen.choose(0, o.dimension - 1)) { (vi : Int) =>
            val i1 = o.linearAssignment(vi, lf).toInterval
            val i2 = o.toInterval.linearAssignment(vi, lf)
            i1 <= i2
          }
        }
      }
    }
  }

  property ("Check that forall X, Y : AbstractOctagon, (X widening Y) >= X, Y (condition 1 of 2 for soundness of widening)") {
    forAll(GenSmallEvenInt) { (d: Int) =>
      forAll(GenMatrix(d)) { (m : M[NonClosed, A]) =>
        val dbmx = e.strongClosure(m)
        val x = new AbstractOctagon[DOM, M, A, B](dbmx, oct, box, e)
        forAll(GenMatrix(d)) { (mm : M[NonClosed, A]) =>
          val dbmy = e.strongClosure(mm)
          val y = new AbstractOctagon[DOM, M, A, B](dbmy, oct, box, e)
          (x widening y) >= x
          (x widening y) >= y
        }
      }
    }
  }

  property ("Check that A U B >= A, B") {
    forAll(GenSmallEvenInt) { (d: Int) =>
      forAll(GenOctagon(d)) { (x : AbstractOctagon[DOM, M, A, B]) =>
        forAll(GenOctagon(d)) { (y : AbstractOctagon[DOM, M, A, B]) =>
          val u = x.union(y)
          u >= x && u >= y
        }
      }
    }
  }

  property ("Check that A \\wedge B <= A, B") {
    forAll(GenSmallEvenInt) { (d: Int) =>
      forAll(GenOctagon(d)) { (x : AbstractOctagon[DOM, M, A, B]) =>
        forAll(GenOctagon(d)) { (y : AbstractOctagon[DOM, M, A, B]) =>
          val i = x.intersection(y)
          i <= x && i <= y
        }
      }
    }
  }


  property ("Check that incrementalClosure is coherent") {
    forAll(GenSmallEvenInt) { (dd: Int) =>
      val d = (dd.abs + 1)*2
      forAll(GenOctagon(d)) { (o : AbstractOctagon[DOM, M, A, B]) =>
        o.linearAssignment(0, DenseLinearForm(Seq(Rational(42))))
        coherent(d, o.dbm)
      }
    }
  }

  property ("Check that incrementalClosure is closed") {
    import utils.singleConstantExactAssignment
    forAll(GenSmallEvenInt) { (dd: Int) =>
      val d = (dd.abs + 1)*2
      forAll(GenOctagon(d)) { (o : AbstractOctagon[DOM, M, A, B]) =>
        o.linearAssignment(0, DenseLinearForm(Seq(Rational(42))))
        closed(d, o.dbm)
      }
    }
  }

  property ("Check that for incrementalClosure m_ij <= (m_{i, bari}+m_{barj, j})/2 holds (condition 3 of 3 for strong closure)") {
    forAll(GenSmallEvenInt) { (dd: Int) =>
      val d = (dd.abs + 1)*2
      forAll(GenOctagon(d)) { (o : AbstractOctagon[DOM, M, A, B]) =>
        o.linearAssignment(0, DenseLinearForm(Seq(Rational(42))))
        condition3(d, o.dbm)
      }
    }
  }

}


class AbstractUtils[A, B <: BoxGenericDomain[A]](
  val box: B, val ifield: InfField[A], GenFinite: Gen[A]) {

  def GenSmallInt : Gen[Int] = Gen.choose(1, 5)
  def GenSmallEvenInt : Gen[Int] = for (n <- GenSmallInt) yield (n * 2)

  val r = new RationalAlgebra()
  def GenRational : Gen[Rational] = for {
    a <- Gen.choose(Int.MinValue, Int.MaxValue)
    b <- Gen.choose(Int.MinValue, Int.MaxValue)
  } yield Rational(a) / Rational(b)

  def GenArbitraryLf(n: Int): Gen[LinearForm] = for {
    coeffs <- Gen.containerOfN[List, Rational](n, GenRational)
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

  // def GenSubsetOf(min: Int, max: Int): Gen[List[Int]] =
  //   if (min > max) Gen.const(Nil) else for {
  //     b <- Gen.oneOf(0, 1)
  //     xs <- GenSubsetOf(min + 1, max)
  //   } yield if (b == 1) min :: xs else xs

  // def GenPartitionOf(list: List[Int]): Gen[List[List[Int]]] = list match {
  //   case Nil => Gen.const(Nil :: Nil)
  //   case x :: Nil => Gen.const((x :: Nil) :: Nil)
  //   case x :: xs => for {
  //     pss <- GenPartitionOf(xs)
  //     b <- Gen.oneOf(0,1)
  //   } yield pss match {
  //     case (p :: ps) => if (b == 1) ((x :: p) :: ps) else ((x :: Nil) :: p :: ps)
  //   }
  // }

  // // // These are for disabling shrinking
  // // // From https://gist.github.com/davidallsopp/f65d73fea8b5e5165fc3
  // // //
  // // // TODO Find a better way

  // // import org.scalacheck.Shrink
  // // implicit val noShrink: Shrink[Int] = Shrink.shrinkAny
  // // implicit val noShrink2: Shrink[(RationalExt, RationalExt)] = Shrink.shrinkAny
  // // implicit val noShrink3: Shrink[FunMatrix[RationalExt]] = Shrink.shrinkAny
  // def HUGE = 1024
  // def SMALL = 1024

  def GenOrderedPair : Gen[(A, A)] =
    for { a <- GenFinite; b <- GenFinite }
    yield
      if (ifield.compare(a,b) == LT || ifield.compare(a,b) == EQ)
        (a, b) else (b, a)

  def GenOrderedDistinctFinitePair : Gen[(A, A)] = {
    val g: Gen[(A, A)] = for {
      low <- GenFinite
      high <- GenFinite
    } yield (low, high)
    g suchThat ({ case (i, j) => ifield.compare(i, j) == LT })
  }

  val GenOrderedDistinctPair =
    GenOrderedPair.suchThat(
      (pair:(A, A)) => ifield.compare(pair._2, pair._1) == GT)

}
