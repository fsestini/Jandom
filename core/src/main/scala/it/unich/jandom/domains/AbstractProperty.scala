/**
 * Copyright 2013 Gianluca Amato <gamato@unich.it>
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

package it.unich.jandom.domains

/**
 * The base class for all abstract properties, i.e. elements of abstract domains. Abstract
 * properties implements a set of poset-like operations, such as union/intersection
 * (corresponding to meet/join) and widening/narrowing. Moreover, abstract properties are
 * partially ordered. Abstract properties use F-bounded polymorhpism to ensure type safety,
 * hence a concrete class `C` implementing an abstract property should inherit from
 * `AbstractProperty[C]`.
 *
 * Properties are partitioned in fibers. Binary operators are guaranteed to work when
 * both elements are part of the same fiber. Finally, properties are immutable.
 *
 * Analogously to `isTop`, `isBottom` and `isEmpty`, comparison between abstract elements
 * should return `true` ONLY IF the relationship actually holds between elements. However,
 * it is not required to return `true` every time the relation holds. The `tryCompare`
 * method, in particular, should return `None` each time it cannot (or does not want)
 * determine the relationship between two abstract properties.
 *
 * @tparam Property the real class we are endowing with the AbstractProperty quality.
 * @define NOTEFIBER `this` and `that` should generally be part of the same fiber.
 * @author Gianluca Amato <gamato@unich.it>
 */
trait AbstractProperty[Property <: AbstractProperty[Property]] extends PartiallyOrdered[Property] {
  /**
   * The class of abstract domains which containts this properties.
   */
  type Domain <: AbstractDomain

  /**
   * Returns the abstract domain corresponding to this property.
   */
  def domain: Domain

  /**
   * The standard widening for two abstract properties.
   * @param that the abstract object to be widened with `this`. `that` is NOT assumed to be bigger than `this`.
   * @note $NOTEFIBER
   * @return the widening of the two abstract properties.
   */
  def widening(that: Property): Property

  /**
   * The standard widening for two abstract properties.
   * @param that the abstract object to be narrowed with `this`. `that` IS assumed to be smaller than `this`.
   * @note $NOTEFIBER
   * @return the narrowing of the two abstract properties.
   */
  def narrowing(that: Property): Property

  /**
   * Compute an upper bound of two abstract properties. If it is possible and convenient, this should compute
   * the least upper bound, but it is not a requirement.
   * @param that the abstract object to join with `this`.
   * @note $NOTEFIBER
   * @return an upper bound of the two abstract properties.
   */
  def union(that: Property): Property

  /**
   * Compute a lower bound of two abstract properties. If it is possible and convenient, this should compute
   * the greatest lower bound, but it is not a requirement.
   * @param that the abstract object to meet with `this`.
   * @note $NOTEFIBER
   * @return a lower bound of the two abstract properties.
   */
  def intersection(that: Property): Property

  /**
   * Returns true ONLY IF this is the top element on the fiber. A top element
   * is bigger than all the other elements, is neutral for intersection and
   * narrowing, and is absorbing for widening and union.
   */
  def isTop: Boolean

  /**
   * Returns true ONLY IF this is the bottom element of the fiber. The
   * opposite is not always true. A bottom element is smaller than all the other elements,
   * is neutral for union and widening and is absorbing for intersection and narrowing.
   */
  def isBottom: Boolean

  /**
   * Returns true ONLY IF this an empty element, i.e. it represents unreachability. If
   * `x.isEmpty` is true, the same happens for `x.isBottom`, but the opposite does
   * not always hold.
   */
  def isEmpty: Boolean

  /**
   * Returns the top property on the same fiber as `this`.
   */
  def top: Property

  /**
   * Returns the bottom property on the same fiber as `this`.
   */
  def bottom: Property

  /**
   * Determines if two properties are the same on the base of `tryCompareTo`
   */
  override def equals(other: Any): Boolean = other match {
    case other: AbstractProperty[Property] => this.tryCompareTo(other) == Some(0)
    case _ => false
  }
}