package it.unich.jandom.domains.numerical
import it.unich.jandom.domains.WideningDescription

class OctagonDomain private[numerical] () extends NumericalDomain {

  // Copied from BoxDoubleDomain.
  val widenings = Seq(WideningDescription.default[Property])

  def apply(): Property = ???


  /**
    * @inheritdoc
    * @note @inheritdoc
    * @throws $ILLEGAL
    */
  def top(dimension: Int): Property = ???

  /**
    * @inheritdoc
    * @note @inheritdoc
    * @throws $ILLEGAL
    */
  def bottom(n: Int): Property = ???

  def mkString(vars: Seq[String]): String = ???

  final class Property extends NumericalProperty[Property] {
    type Domain = OctagonDomain

    def domain = OctagonDomain.this

    private def normalized: Boolean = ???

    /**
      * @inheritdoc
      * @note @inheritdoc
      * @throws $ILLEGAL
      */
    def widening(that: Property): Property = ???

    /**
      * @inheritdoc
      * @note @inheritdoc
      * @throws $ILLEGAL
      */
    def narrowing(that: Property): Property = ???

    /**
      * @inheritdoc
      * It is equivalent to `intersectionWeak`.
      * @note @inheritdoc
      * @throws $ILLEGAL
      */
    def intersection(that: Property): Property = intersectionWeak(that)

    /**
      * @inheritdoc
      * The union of two parallelotopes is not a parallelotope. Moreover, there is no least
      * parallelotopes containing the union of two parallelotopes. Hence, this methods uses
      * heuristics to find a good result.
      * @note @inheritdoc
      * @throws $ILLEGAL
      */
    def union(that: Property): Property = ???

    /**
      * This is a variant of `union` using weak join. The shape of the resulting
      * parallelotope is the same shap of `this`.
      *
      * @param that the abstract object to be joined with `this`.
      * @note $NOTEDIMENSION
      * @return the weak union of the two abstract objects.
      */
    def unionWeak(that: Property): Property = ???

    /**
      * This is the weak intersection of two abstract objects. The shape of the resulting
      * parallelotope is the same shap of `this`.
      *
      * @param that the abstract object to be intersected with `this`.
      * @note $NOTEDIMENSION
      * @return the intersection of the two abstract objects.
      */
    def intersectionWeak(that: Property): Property = ???

    /**
      * @inheritdoc
      * @note @inheritdoc
      * @todo @inheritdoc
      * @throws $ILLEGAL
      */
    def linearAssignment(n: Int, lf: LinearForm): Property = ???

    /**
      * @inheritdoc
      * @note @inheritdoc
      * @todo @inheritdoc
      * @throws $ILLEGAL
      */
    def linearInequality(lf: LinearForm): Property = ???

    /**
      * @inheritdoc
      * @note @inheritdoc
      * @throws $ILLEGAL
      */
    def linearDisequality(lf: LinearForm): Property = ???

    /**
      * @inheritdoc
      * @note @inheritdoc
      * @throws $ILLEGAL
      */
    def nonDeterministicAssignment(n: Int): Property = ???

    def addVariable(): Property = ???

    def constraints = ???

    def isPolyhedral = true

    /**
      * @inheritdoc
      * @note @inheritdoc
      * @throws $ILLEGAL
      */
    def delVariable(n: Int): Property = ???

    /**
      * @inheritdoc
      */
    def mapVariables(rho: Seq[Int]): Property = ???

    /**
      * Compute the minimum and maximum value of a linear form in a parallelotope.
      *
      * @todo should be generalized to linear forms over arbitrary types.
      * @return a tuple with two components: the first component is the least value, the second component is the greatest value
      * of the linear form over the box.
      */
    //def linearEvaluation(lf: LinearForm): (RationalExt, RationalExt) =

    def minimize(lf: LinearForm) = ???

    def maximize(lf: LinearForm) = ???

    def frequency(lf: LinearForm) = ???

    def dimension = ???.asInstanceOf[Int]

    def isEmpty = ???

    def isTop = ???

    def isBottom = ???

    def bottom = OctagonDomain.this.bottom(dimension)

    def top = OctagonDomain.this.top(dimension)

    def tryCompareTo[B >: Property](that: B)(implicit arg0: (B) => PartiallyOrdered[B]): Option[Int] = ???


    def <=[B >: Property](that: Property)(implicit arg0: (B) => PartiallyOrdered[B]): Boolean = ???

    def >=[B >: Property](that: Property)(implicit arg0: (B) => PartiallyOrdered[B]): Boolean =
      that <= this

    def <[B >: Property](that: Property)(implicit arg0: (B) => PartiallyOrdered[B]): Boolean =
      (this <= that) && !(this >= that)

    def >[B >: Property](that: Property)(implicit arg0: (B) => PartiallyOrdered[B]): Boolean =
      (this >= that) && !(this <= that)

    def mkString(vars: Seq[String]): String = ???
  }

}