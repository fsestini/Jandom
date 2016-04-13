/**
 * Copyright 2013, 2016 Jandom Team
 *
 * This file is part of JANDOM: JVM-based Analyzer for Numerical DOMains
 * JANDOM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * JANDOM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty ofa
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with JANDOM.  If not, see <http://www.gnu.org/licenses/>.
 */

package it.unich.jandom.domains.numerical

import it.unich.jandom.utils.numberext.RationalExt
import it.unich.jandom.domains.CachedTopBottom

/**
 * This is the domain of boxes, also known as the interval domain. Bounds are represented by rational numbers.
 *
 * This domain is implemented natively in Scala and is a sort of reference implementation for numerical domains.
 * Real domains should be better implemented within $PPL or $APRON.
 * @author Gianluca Amato <gianluca.amato@unich.it>
 * @author Marco Rubino <marco.rubino@unich.it>
 *
 * @constructor Builds a box domain uses rational as bounds
 */

class BoxRationalDomain private extends NumericalDomain {
  /**
   * This is the class representing a single box.
   *
   * @constructor Creates a box with given lower and upper bounds.
   * @param low the lower bounds of the box.
   * @param high the upper bounds of the box.
   * @param isEmpty is true when the box is empty. It is needed for the case of 0-dimensional boxes.
   * @note `low`, `high` and `isEmpty` should be normalized according to the method `normalized`
   * @throws IllegalArgumentException if parameters are not correct.
   */

  final class Property(val low: Array[RationalExt], val high: Array[RationalExt], val isEmpty: Boolean) extends NumericalProperty[Property] {
    require(normalized, s"The parameters low: ${low.mkString(",")}, high: ${high.mkString(",")} and isEmpty: ${isEmpty} are not normalized")

    type Domain = BoxRationalDomain

    def domain = BoxRationalDomain.this

    /**
     * This checks whether the box is normalized. This should always be the case. A box is normalized when
     * the lower and higher bounds are of the same length, and either
     *   1. there are no lower bounds equal to +Inf, there are no upper bounds equal to -Inf,
     *       the lower bounds are smaller of the corresponding upper bounds, isEmpty is false, or
     *   2. all the lower bound are +Inf and all the upper bounds are -Inf, isEmpty is true.
     * @return whether the box is normalized.
     */
    private def normalized: Boolean =
      low.length == high.length &&
        (
          low.forall { (x) => !(x.isPosInfinity) } &&
          high.forall { (x) => !(x.isNegInfinity) } &&
          (low, high).zipped.forall(_ <= _) &&
          !isEmpty
          ||
          low.forall { _.isPosInfinity } &&
          high.forall { _.isNegInfinity } &&
          isEmpty)

    /**
     * Return the dot product of `x` and `y`. If element `x(i)` is zero, then `x(i)*y(i)` is
     * `0` independently from the value of `y(i)`.
     * If `remove` is a valid index in `x` and `y`, the factor `x(remove) * y(remove)` is
     * removed from the dot product.
     */
    private def dotprod(x: Seq[RationalExt], y: Seq[RationalExt], remove: Int = -1): RationalExt = {
      var sum: RationalExt = RationalExt.zero
      for (i <- x.indices if i != remove if x(i) != RationalExt.zero) sum = (sum + (x(i) * y(i)))
      sum
    }

    /**
     * @inheritdoc
     * @note @inheritdoc
     * @throws $ILLEGAL
     */
    def union(that: Property): Property = {
      require(dimension == that.dimension)
      val newlow = (this.low, that.low).zipped.map(_ min _)
      val newhigh = (this.high, that.high).zipped.map(_ max _)
      new Property(newlow, newhigh, isEmpty && that.isEmpty)
    }

    /**
     * @inheritdoc
     * @note @inheritdoc
     * @throws $ILLEGAL
     */
    def intersection(that: Property): Property = {
      require(dimension == that.dimension)
      val newlow = (this.low, that.low).zipped.map(_ max _)
      val newhigh = (this.high, that.high).zipped.map(_ min _)
      BoxRationalDomain.this(newlow, newhigh)
    }

    /**
     * This is the standard widening on boxes based on [[http://www.di.ens.fr/~cousot/COUSOTpapers/ISOP76.shtml CC76]].
     * @note @inheritdoc
     * @throws $ILLEGAL
     */
    def widening(that: Property) = {
      require(dimension == that.dimension)
      val newlow = (low, that.low).zipped.map((l1, l2) => if (l1 == RationalExt.PositiveInfinity) l2 else if (l1 <= l2) l1 else RationalExt.NegativeInfinity)
      val newhigh = (high, that.high).zipped.map((l1, l2) => if (l1 == RationalExt.NegativeInfinity) l2 else if (l1 >= l2) l1 else RationalExt.PositiveInfinity)
      new Property(newlow, newhigh, isEmpty && that.isEmpty)
    }

    /**
     * This is the standard narrowing on boxes based on [[http://www.di.ens.fr/~cousot/COUSOTpapers/ISOP76.shtml CC76]].
     * @note @inheritdoc
     * @throws $ILLEGAL
     */
    def narrowing(that: Property) = {
      require(dimension == that.dimension)

      if (that.isEmpty) {
        that
      } else {
        val newlow = (low, that.low).zipped.map((l1, l2) => if (l1 == RationalExt.NegativeInfinity) l2 else l1 min l2)
        val newhigh = (high, that.high).zipped.map((l1, l2) => if (l1 == RationalExt.PositiveInfinity) l2 else l1 max l2)
        BoxRationalDomain.this(newlow, newhigh)
      }
    }

    /**
     * Compute the minimum and maximum value of a linear form in a box.
     * @param lf a linear form.
     * @return a tuple with two components: the first component is the least value, the second component is the greatest value
     * of the linear form over the box.
     */
    private def linearEvaluation_r(lf: LinearForm[RationalExt]): (RationalExt, RationalExt) = {
      require(lf.dimension <= dimension)
      var newlow: RationalExt = lf.known
      var newhigh: RationalExt = lf.known
      val coeffs = lf.homcoeffs
      if (isEmpty && coeffs.exists { _ != RationalExt.zero })
        (RationalExt.PositiveInfinity, RationalExt.NegativeInfinity)
      else {
        for (i <- coeffs.indices) {
          if (coeffs(i) < RationalExt.zero) {
            newlow = (newlow + (coeffs(i) * high(i)))
            newhigh = (newhigh + (coeffs(i) * low(i)))
          } else if (coeffs(i) > RationalExt.zero) {
            newlow = (newlow + (coeffs(i) * low(i)))
            newhigh = (newhigh + (coeffs(i) * high(i)))
          }
        }
        (newlow, newhigh)
      }
    }

    def minimize_r(lf: LinearForm[RationalExt]) = linearEvaluation_r(lf)._1

    def maximize_r(lf: LinearForm[RationalExt]) = linearEvaluation_r(lf)._2

    def frequency_r(lf: LinearForm[RationalExt]) = {
      val (min, max) = linearEvaluation_r(lf)
      if (min == max) Some(min) else None
    }

    // The following three methods are only needed because the Jandom API
    // wants linear forms over doubles.
    def linearEvaluation(lf: LinearForm[Double]): (Double, Double) = {
      val le = linearEvaluation_r(LinearForm(lf.coeffs map { RationalExt(_) }: _*))
      (le._1.toDouble, le._2.toDouble)
    }

    def minimize(lf: LinearForm[Double]) = linearEvaluation(lf)._1

    def maximize(lf: LinearForm[Double]) = linearEvaluation(lf)._2

    def frequency(lf: LinearForm[Double]) = {
      val (min, max) = linearEvaluation(lf)
      if (min == max) Some(min) else None
    }

    /**
     * Compute the corner of the box which minimizes a linear form.
     * todo should be generalized to linear forms over arbitrary types.
     * @param coeff the homogeneous coefficients.
     * @return the coordinates of the point which minimizes the linear form.
     */
    private def linearArgmin(lf: LinearForm[RationalExt]): Seq[RationalExt] = {
      require(lf.dimension <= dimension)
      (lf.homcoeffs.zipWithIndex) map { case (c, i) => if (c > RationalExt.zero) low(i) else high(i) }
    }

    /**
     * Compute the corner of the box which maximizes a linear form.
     * @todo should be generalized to linear forms over arbitrary types.
     * @param coeff the homogeneous coefficients
     * @return the coordinates of the point which maximizes the linear form
     */
    private def linearArgmax(lf: LinearForm[RationalExt]): Seq[RationalExt] = {
      require(lf.dimension <= dimension)
      (lf.homcoeffs.zipWithIndex) map { case (c, i) => if (c < RationalExt.zero) low(i) else high(i) }
    }

    /**
     * @inheritdoc
     * @note @inheritdoc
     * @throws $ILLEGAL
     */
    def nonDeterministicAssignment(n: Int): Property = {
      require(n < low.length && n >= 0)
      if (isEmpty)
        this
      else
        new Property(low.updated(n, RationalExt.NegativeInfinity), high.updated(n, RationalExt.PositiveInfinity), false)
    }

    def linearAssignment(n: Int, lf: LinearForm[Double]): Property = linearAssignment_m(n, LinearForm(lf.coeffs map { RationalExt(_) }: _*))

    /**
     * @inheritdoc
     * @note @inheritdoc
     * @todo @inheritdoc
     * @throws $ILLEGAL
     */
    def linearAssignment_m(n: Int, lf: LinearForm[RationalExt]): Property = {
      require(n < low.length && n >= 0 && lf.dimension <= dimension)
      if (isEmpty)
        this
      else {
        val interval = linearEvaluation_r(lf)
        new Property(low.updated(n, interval._1), high.updated(n, interval._2), false)
      }

    }

    /**
     * @inheritdoc
     * @note @inheritdoc
     * @todo @inheritdoc
     * @throws $ILLEGAL
     */
    def linearInequality(lf: LinearForm[Double]): Property =
      linearInequality_r(LinearForm(lf.coeffs map { RationalExt(_) }: _*))

    def linearInequality_r(lf: LinearForm[RationalExt]): Property = {
      require(lf.dimension <= dimension)

      /* if the box is empty the result is empty */
      if (isEmpty) return this

      /* check if result is empty */
      val lfArgmin = linearArgmin(lf)
      val lfMin = linearEvaluation_r(lf)._1
      if (lfMin > RationalExt.zero) return BoxRationalDomain.this.bottom(dimension)

      val newlow = low.clone
      val newhigh = high.clone

      val coeffs = lf.homcoeffs
      val known = lf.known

      val infinities = (coeffs.indices) filter { i => lfArgmin(i).isInfinity && coeffs(i) != RationalExt.zero }
      infinities.size match {
        case 0 =>
          for (i <- coeffs.indices) {
            if (coeffs(i) < RationalExt.zero) newlow(i) = low(i) max (lfArgmin(i) - lfMin / coeffs(i))
            if (coeffs(i) > RationalExt.zero) newhigh(i) = high(i) min (lfArgmin(i) - lfMin / coeffs(i))
          }
        case 1 => {
          val posinf = infinities.head
          if (coeffs(posinf) < RationalExt.zero)
            newlow(posinf) = low(posinf) max ((-dotprod(coeffs, lfArgmin, posinf) - known) / coeffs(posinf))
          else
            newhigh(posinf) = high(posinf) min ((-dotprod(coeffs, lfArgmin, posinf) - known) / coeffs(posinf))
        }
        case _ =>
      }
      BoxRationalDomain.this(newlow, newhigh)
    }

    def constraints = {
      if (isEmpty)
        Seq(LinearForm(1))
      else {
        val set1 = for (i <- 0 until dimension; if !low(i).isInfinity) yield -LinearForm.v[RationalExt](i) + low(i)
        val set2 = for (i <- 0 until dimension; if !high(i).isInfinity) yield LinearForm.v[RationalExt](i) - high(i)
        (set1 ++ set2) map { _.toDouble }
      }
    }

    def isPolyhedral = true

    /**
     * @inheritdoc
     * @note @inheritdoc
     * @throws $ILLEGAL
     */
    def linearDisequality(lf: LinearForm[Double]): Property =
      linearDisequality_r(LinearForm(lf.coeffs map { RationalExt(_) }: _*))

    def linearDisequality_r(lf: LinearForm[RationalExt]): Property = {
      val count = lf.homcoeffs.count(_ != RationalExt.zero)
      count match {
        case 0 =>
          if (lf.known == RationalExt.zero) bottom else this
        case 1 =>
          val dim = lf.homcoeffs.indexWhere(_ != RationalExt.zero)
          if (low(dim) == lf.known && high(dim) == lf.known)
            bottom
          else
            this
        case _ => this
      }
    }

    def addVariable: Property =
      if (isEmpty)
        BoxRationalDomain.this.bottom(dimension + 1)
      else
        BoxRationalDomain.this(low :+ RationalExt.NegativeInfinity, high :+ RationalExt.PositiveInfinity)

    /**
     * @inheritdoc
     * This is a complete operator for boxes.
     * @note @inheritdoc
     * @throws $ILLEGAL
     */
    def delVariable(n: Int): Property = {
      require(n < low.length && n >= 0)
      val newlow = new Array[RationalExt](dimension - 1)
      val newhigh = new Array[RationalExt](dimension - 1)
      Array.copy(low, 0, newlow, 0, n)
      Array.copy(high, 0, newhigh, 0, n)
      Array.copy(low, n + 1, newlow, n, dimension - n - 1)
      Array.copy(high, n + 1, newhigh, n, dimension - n - 1)
      new Property(newlow, newhigh, isEmpty)
    }

    /**
     * @inheritdoc
     * This is a complete operator for boxes.
     * @note @inheritdoc
     * @throws IllegalArgumentException if parameters are not correct (but we do not check injectivity of `rho`)
     */
    def mapVariables(rho: Seq[Int]) = {
      require(rho.length == dimension)
      val newdim = rho.count(_ >= 0)
      require(rho forall { i => i >= -1 && i < newdim })
      // we do not check injectivity
      val newlow = new Array[RationalExt](newdim)
      val newhigh = new Array[RationalExt](newdim)
      for ((newi, i) <- rho.zipWithIndex; if newi >= 0) {
        newlow(newi) = low(i)
        newhigh(newi) = high(i)
      }
      new Property(newlow, newhigh, isEmpty)
    }

    /**
     * @inheritdoc
     * @throws $ILLEGAL
     */
    def mkString(vars: Seq[String]): String = {
      require(vars.length >= dimension)
      if (isEmpty)
        "empty"
      else {
        val bounds = for (i <- 0 until dimension) yield {
          if (low(i) < high(i))
            low(i) + " <= " + vars(i) + " <= " + high(i)
          else vars(i) + " = " + high(i)
        }
        bounds.mkString("[ ", " , ", " ]")
      }
    }

    val dimension: Int = low.length

    def isBottom = isEmpty

    def isTop = !isEmpty && low.forall(_.isNegInfinity) && high.forall(_.isPosInfinity)

    def bottom = BoxRationalDomain.this.bottom(low.length)

    def top = BoxRationalDomain.this.top(low.length)

    def tryCompareTo[B >: Property](other: B)(implicit arg0: (B) => PartiallyOrdered[B]): Option[Int] = other match {
      case other: Property =>
        require(dimension == other.dimension)
        (isEmpty, other.isEmpty) match {
          case (true, true) => Some(0)
          case (false, true) => Some(1)
          case (true, false) => Some(-1)
          case (false, false) =>
            val lowpairs = (this.low, other.low).zipped
            val highpairs = (this.high, other.high).zipped
            if (lowpairs.forall(_ == _) && highpairs.forall(_ == _))
              Some(0)
            else if (lowpairs.forall(_ <= _) && highpairs.forall(_ >= _))
              Some(1)
            else if (lowpairs.forall(_ >= _) && highpairs.forall(_ <= _))
              Some(-1)
            else
              None
        }
      case _ => None
    }

    // override def hashCode: Int = 41 * (41 + low.hashCode) + high.hashCode
  }

  /**
   * Returns a normalized box with given bounds.
   * @param low lower bounds.
   * @param high upper bounds.
   * @note `low` should have the same length as `high`.
   * @return the normalized box with the specified bounds.
   * @throws $ILLEGAL
   */
  def apply(low: Array[RationalExt], high: Array[RationalExt]): Property = {
    require(low.length == high.length)
    if ((low, high).zipped.exists(_ > _))
      bottom(low.length)
    else
      new Property(low, high, false)
  }

  /**
   * Returns a box consisting of the single point `poDouble`.
   */
  def apply(poDouble: Array[RationalExt]): Property = apply(poDouble, poDouble)

  /**
   * @inheritdoc
   * @note @inheritdoc
   * @throws $ILLEGAL
   */
  def top(n: Int): Property =
    new Property(Array.fill(n)(RationalExt.NegativeInfinity), Array.fill(n)(RationalExt.PositiveInfinity), false)

  /**
   * @inheritdoc
   * @note @inheritdoc
   * @throws $ILLEGAL
   */
  def bottom(n: Int): Property =
    new Property(Array.fill(n)(RationalExt.PositiveInfinity), Array.fill(n)(RationalExt.NegativeInfinity), true)
}

object BoxRationalDomain {
  /**
   * Returns an abstract domain for boxes with rational bounds.
   */
  def apply() = domain

  /**
   * The domain of boxes with rational bounds
   */
  private val domain = new BoxRationalDomain with CachedTopBottom
}
