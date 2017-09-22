package it.unich.jandom.domains.numerical.octagon
import scalaz._

import it.unich.jandom.domains.numerical.octagon.variables.Dimension
import it.unich.jandom.domains.numerical.octagon.variables.CountOps._

import scala.language.higherKinds

// Simple implementation of a square matrix by means of nested vectors.
case class VecMatrix[A](private val vec: Vector[Vector[A]], val dimension: Dimension) {

  def update(i: Int, j: Int, x: A): VecMatrix[A] = {
    require(inDimension(i, j, dimension))
    new VecMatrix(vec.updated(i, vec(i).updated(j, x)), dimension)
  }

  def update(updater: (Int, Int) => A): VecMatrix[A] = {
    val range = allIndices(dimension).toVector
    val mat = range.map(i => range.map(j => updater(i, j)))
    new VecMatrix(mat, dimension)
  }

  def apply(i: Int, j: Int): A = {
    require(inDimension(i, j, dimension))
    vec(i)(j)
  }

  def combine[B, C](that: VecMatrix[B], f: (A, B) => C): VecMatrix[C] = {
    require(this.dimension == that.dimension)
    val mat = (this.vec zip that.vec) map { case (row1, row2) =>
      (row1 zip row2) map { case (elem1, elem2) =>
        f(elem1, elem2)
      }
    }

    new VecMatrix(mat, dimension)
  }

  def toList: List[A] = vec.flatten.toList
}

object VecMatrix {

  def apply[A](dim: Dimension, x: A): VecMatrix[A] =
    new VecMatrix[A](fill(dim)(x), dim)

  def apply[A](dim: Dimension, f: (Int, Int) => A): VecMatrix[A] = {
    new VecMatrix(
      (allIndices(dim)).map(
        (i: Int) => (allIndices(dim)).map(
          (j: Int) => f(i,j)).toVector).toVector, dim)
  }

}

object VecMatrixMatrixInstance {
  implicit val vecMatrixIsMatrix: Matrix[VecMatrix] = new Matrix[VecMatrix] {
    def combine[A, B, C](f: (A, B) => C)
                        (ma: VecMatrix[A], mb: VecMatrix[B]): VecMatrix[C] =
      ma.combine(mb, f)
    def update[A](i: Int, j: Int, x: A)(m: VecMatrix[A]): VecMatrix[A] =
      m.update(i, j, x)
    def update[A](updater: (Int, Int) => A)(m: VecMatrix[A]): VecMatrix[A] =
      m.update(updater)
    def get[A](i: Int, j: Int)(m: VecMatrix[A]): A = m(i, j)
    def foldMap[A, B](fa: VecMatrix[A])(f: (A) => B)(implicit F: Monoid[B]): B =
      fa.toList.map(f).fold(F.zero)((x, y) => F.append(x, y))

    def foldRight[A, B](fa: VecMatrix[A], z: => B)(f: (A, => B) => B): B =
      fa.toList.foldRight(z)((x, y) => f(x, y))

    def pure[A](dim: Dimension, x: A): VecMatrix[A] = VecMatrix(dim, x)
    def make[A](f: ((Int, Int) => A), dim: Dimension) = VecMatrix(dim, f)

    def dimension[A](m: VecMatrix[A]): Dimension = m.dimension

  }
}
