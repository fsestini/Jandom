package it.unich.jandom.domains.numerical.octagon.testutils
import it.unich.jandom.domains.numerical.octagon._
import it.unich.jandom.domains.numerical._
import it.unich.jandom.domains.numerical.octagon.fast._

import variables._

import org.scalacheck.Gen
import org.scalacheck.Arbitrary
import org.scalacheck.Arbitrary.arbitrary
import spire.math.Rational
import spire.math.RationalAlgebra
import it.unich.jandom.utils.numberext.RationalExt

class Utils(val box: BoxRationalDomain) {
  val me =  FunMatrixMatrixInstance.funMatrixIsMatrix
  type FunDBM[B,C] = (DBM[FunMatrix, B, C])
  // val box = BoxRationalDomain()

  val r = new RationalAlgebra()

  def GenSmallInt : Gen[Int] = Gen.choose(1, 10)
  def GenSmallEvenInt : Gen[Int] = for (n <- GenSmallInt) yield (n * 2)
  def GenInf : Gen[RationalExt] = Gen.const(RationalExt.PositiveInfinity)
  def GenRational : Gen[Rational] = for {
    a <- Gen.choose(Int.MinValue, Int.MaxValue)
    b <- Gen.choose(Int.MinValue, Int.MaxValue)
  } yield        Rational(a) / Rational(b)

  def GenFiniteRationalExt = for (a <- GenRational) yield RationalExt(a)
  def GenFiniteRationalExtsAndInf(pInf: Int) : Gen[RationalExt] = Gen.frequency(
    (pInf, GenInf),
    (100 - pInf, GenFiniteRationalExt)
  )

  def GenArbitraryLf(n: Int): Gen[LinearForm] = for
  {
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

  // // These are for disabling shrinking
  // // From https://gist.github.com/davidallsopp/f65d73fea8b5e5165fc3
  // //
  // // TODO Find a better way

  // import org.scalacheck.Shrink
  // implicit val noShrink: Shrink[Int] = Shrink.shrinkAny
  // implicit val noShrink2: Shrink[(RationalExt, RationalExt)] = Shrink.shrinkAny
  // implicit val noShrink3: Shrink[FunMatrix[RationalExt]] = Shrink.shrinkAny
  def HUGE = 1024
  def SMALL = 1024

  def GenOrderedPair : Gen[(RationalExt, RationalExt)] = for {
    a <- GenFiniteRationalExt
    b <- GenFiniteRationalExt
  } yield if (a <= b)
    (a, b)
  else (b, a)


  def GenOrderedFinitePair : Gen[(Rational, Rational)] = for {
    low <- Gen.choose(Int.MinValue, Int.MaxValue)
    high <- Gen.choose(low, Int.MaxValue)
  } yield (low, high)

  val GenOrderedDistinctPair = GenOrderedPair.suchThat((pair:(RationalExt, RationalExt)) =>pair._2 > pair._1)
  val GenOrderedDistinctFinitePair = GenOrderedFinitePair.suchThat((pair:(Rational, Rational)) =>pair._2 > pair._1)

  def GenFunMatrix(d: Int, pTop: Int = 20, pInf: Int = 20) : Gen[FunMatrix[RationalExt]] = for {
    rowSeq <- Gen.containerOfN[Array, RationalExt](d, GenFiniteRationalExtsAndInf(pInf))
    arrayOfRows <- Gen.containerOfN[Array, Array[RationalExt]](d, rowSeq)
  } yield new FunMatrix[RationalExt](
      ((i: Int, j: Int) =>
        if (i == j) 0
        else arrayOfRows(i)(j)), Dimension(d))

  def GenTop(nOfVars: Int): Gen[FunDBM[Closed, RationalExt]] = Gen.const(new ClosedDBM[FunMatrix, RationalExt](
    // Caveat: top is not strongly closed from the defn, but we have it as an instance of FunDBM[Closed,_].
        FunMatrix[RationalExt]((i: Int, j: Int) =>
        if (i == j) 0 // If (i,i) != 0 we have bottom
        else RationalExt.PositiveInfinity, Dimension(nOfVars*2)))(me)
  )

  def GenClosedFunDBMOrTop(nOfVars: Int, pBot: Int = 10, pTop: Int = 20, pInf: Int = 30) : Gen[FunDBM[Closed, RationalExt]] =
    Gen.frequency(
      (100 - pTop, GenTop(nOfVars)),
      (pTop, GenClosedFunDBM(nOfVars, pInf))
    )


  def GenClosedFunDBM(nOfVars: Int, pInf: Int = 30) : Gen[FunDBM[Closed, RationalExt]] =
      for { m <- GenFunMatrix(nOfVars * 2, pInf) }
      yield
      {
        // caveat: there is no guarantee re: the distribution of bottoms, should probably include a few pre-computed non-bottom ones?
        val closure = (new BagnaraStrongClosure[FunMatrix]()(me)).strongClosure(m)
        if (closure == None) {
          new BottomDBM(VarCount(nOfVars))
        } else {
          new ClosedDBM[FunMatrix, RationalExt](closure.get)(me)
        }
      }

  def GenHalfMatrix(d: Int, pInf: Int = 20): Gen[HalfMatrix[RationalExt]] = for {
    rowSeq <- Gen.containerOfN[Array, RationalExt](d, GenFiniteRationalExtsAndInf(pInf))
    arrayOfRows <- Gen.containerOfN[Array, Array[RationalExt]](d, rowSeq)
  } yield
    HalfMatrix[RationalExt](
      (i: Int, j: Int) => if (i == j) 0 else arrayOfRows(i)(j),
      VarCount(d / 2))

}
