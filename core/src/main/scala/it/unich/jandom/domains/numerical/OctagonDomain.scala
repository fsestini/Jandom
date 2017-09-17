package it.unich.jandom.domains.numerical
import it.unich.jandom.domains.WideningDescription
import it.unich.jandom.domains.numerical.octagon._
import it.unich.jandom.domains.numerical.octagon.variables.VarCount
import scala.language.higherKinds
import scala.reflect.ClassTag
import it.unich.jandom.utils.numberext.RationalExt

class OctagonDomain[M[_,N], N, B <: BoxGenericDomain[N]] (
  e: DifferenceBoundMatrix[M] { type PosetConstraint[A] = InfField[A]}, b: B)
  (implicit ifield: InfField[N], ratToN: (RationalExt => N), ct: ClassTag[N])
    extends NumericalDomain {

  // Copied from BoxDoubleDomain.
  val widenings = Seq(WideningDescription.default[Property])

  /**
    * @inheritdoc
    * @note @inheritdoc
    * @throws $ILLEGAL
    */
  def top(dimension: Int): Property = AbstractOctagon[OctagonDomain[M, N, B], M, N, B](e.topDBM(VarCount(dimension)), this, b, e)

  /**
    * @inheritdoc
    * @note @inheritdoc
    * @throws $ILLEGAL
    */
  def bottom(dimension: Int): Property = AbstractOctagon[OctagonDomain[M, N, B], M, N, B](e.bottomDBM(VarCount(dimension)), this, b, e)

  type Property = AbstractOctagon[OctagonDomain[M, N, B], M, N, B]
}

object OctagonDomain {
  val VecDBMInstance = (new DBMInstance[VecMatrix]()(VecMatrixMatrixInstance.vecMatrixIsMatrix))
  def id (x: RationalExt): RationalExt = x
  import InfField.ifieldRationalExt
  type VecDBM[S,A] = DBM[VecMatrix, S, A]
  val b =  BoxRationalDomain()
  def apply() = new OctagonDomain[VecDBM, RationalExt, BoxRationalDomain](VecDBMInstance.funDBM, b)
}
