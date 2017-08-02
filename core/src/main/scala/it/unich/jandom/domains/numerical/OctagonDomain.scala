package it.unich.jandom.domains.numerical
import it.unich.jandom.domains.WideningDescription
import it.unich.jandom.domains.numerical.octagon._
import scala.language.higherKinds

class OctagonDomain[M[_,_]] private[numerical] (
  e: DifferenceBoundMatrix[M] { type PosetConstraint[A] = InfField[A]})
    extends NumericalDomain {

  // Copied from BoxDoubleDomain.
  val widenings = Seq(WideningDescription.default[Property])

  /**
    * @inheritdoc
    * @note @inheritdoc
    * @throws $ILLEGAL
    */
  def top(dimension: Int): Property = AbstractOctagon(e.topDBM(dimension), this, e)

  /**
    * @inheritdoc
    * @note @inheritdoc
    * @throws $ILLEGAL
    */
  def bottom(dimension: Int): Property = AbstractOctagon(e.bottomDBM(dimension), this, e)

  type Property = AbstractOctagon[OctagonDomain[M], M]
}

object OctagonDomain {
  def apply() = new OctagonDomain(FunDBMInstance.funDBM)
}
