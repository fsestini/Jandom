package it.unich.jandom.domains.numerical.octagon

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

class RationalVecMatrixSpecification extends AbstractOctagonSpecification {

  type A = RationalExt
  type B = BoxRationalDomain
  type M[S,X] = DBM[VecMatrix, S, X]

  def GenRational: Gen[Rational] = for {
    a <- Gen.choose(Int.MinValue, Int.MaxValue)
    b <- Gen.choose(Int.MinValue, Int.MaxValue)
  } yield Rational(a) / Rational(b)

  def GenFinite = for (a <- GenRational) yield RationalExt(a)

  def GenFiniteRationalExtsAndInf(pInf: Int) : Gen[RationalExt] =
    Gen.frequency(
      (pInf, Gen.const(RationalExt.PositiveInfinity)),
      (100 - pInf, GenFinite))

  def GenMatrix(d: Int, pInf: Int = 20): Gen[M[NonClosed, RationalExt]] = for {
    rowSeq <- Gen.containerOfN[Array, RationalExt](d, GenFiniteRationalExtsAndInf(pInf))
    arrayOfRows <- Gen.containerOfN[Array, Array[RationalExt]](d, rowSeq)
  } yield {
    val f: (Int, Int) => RationalExt = (i, j) => if (i == j) 0 else arrayOfRows(i)(j)
    e.fromFun(Dimension(d), f)
  }

  import Utils._

  lazy val box = BoxRationalDomain()
  def toRat = (x: RationalExt) => x.value

  def createInterval(low: Array[RationalExt], high: Array[RationalExt]) =
    new box.Property(low, high, false)

  implicit def ct: ClassTag[A] = ClassTag(classOf[RationalExt])
  def e =
    (new DBMInstance[VecMatrix]()(VecMatrixMatrixInstance.vecMatrixIsMatrix)).funDBM

  implicit val ifield: InfField[RationalExt] = InfField.ifieldRationalExt
  implicit val ratToN: RationalExt => RationalExt = (x: RationalExt) => x

}
