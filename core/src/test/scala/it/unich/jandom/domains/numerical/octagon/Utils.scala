package it.unich.jandom.domains.numerical.octagon.testutils
import it.unich.jandom.domains.numerical.octagon._
import it.unich.jandom.domains.numerical._

import variables._

import org.scalacheck.Gen
import org.scalacheck.Arbitrary
import org.scalacheck.Arbitrary.arbitrary
import spire.math.Rational
import spire.math.RationalAlgebra

object Utils {
  val me =  FunMatrixMatrixInstance.funMatrixIsMatrix
  type FunDBM[B,C] = (DBM[FunMatrix, B, C])
  val box = BoxDoubleDomain(false)
  val r = new RationalAlgebra()
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
        pairs : Array[(Double, Double)] <- Gen.containerOfN[Array, (Double, Double)](n, GenOrderedPair)
      } yield (new box.Property(pairs.unzip._1, pairs.unzip._2, false))
    }

  def GenSmallInt : Gen[Int] = Gen.choose(1, 5)
  def GenSmallEvenInt : Gen[Int] = for (n <- GenSmallInt) yield (n * 2)
  def GenInf : Gen[Double] = Gen.const(Double.PositiveInfinity)
  def GenDoublesAndInf(pInf: Int) : Gen[Double] = Gen.frequency(
    (pInf, GenInf),
    (100 - pInf, Gen.choose(Float.MinValue.toDouble, Float.MaxValue.toDouble))
    // We generate only single precision floats to avoid false positives due to 754 edge cases
  )

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

  // // These are for disabling shrinking
  // // From https://gist.github.com/davidallsopp/f65d73fea8b5e5165fc3
  // //
  // // TODO Find a better way

  // import org.scalacheck.Shrink
  // implicit val noShrink: Shrink[Int] = Shrink.shrinkAny
  // implicit val noShrink2: Shrink[(Double, Double)] = Shrink.shrinkAny
  // implicit val noShrink3: Shrink[FunMatrix[Double]] = Shrink.shrinkAny

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

  def GenTop(nOfVars: Int): Gen[FunDBM[Closed, Double]] = Gen.const(new ClosedDBM[FunMatrix, Double](
    // Caveat: top is not strongly closed from the defn, but we have it as an instance of FunDBM[Closed,_].
        FunMatrix[Double]((i: Int, j: Int) =>
        if (i == j) 0 // If (i,i) != 0 we have bottom
        else Double.PositiveInfinity, Dimension(nOfVars*2)))(me)
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
        val closure = (new BagnaraStrongClosure[FunMatrix, Double]()(me)).strongClosure(m)
        if (closure == None) {
          new BottomDBM(VarCount(nOfVars))
        } else {
          new ClosedDBM[FunMatrix, Double](closure.get)(me)
        }
      }

  def checkIsLegal(m : FunMatrix[Double]) : Boolean =
    CountOps.allIndices(m.dimension).forall(
      (i)=> CountOps.allIndices(m.dimension).forall(
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
}
