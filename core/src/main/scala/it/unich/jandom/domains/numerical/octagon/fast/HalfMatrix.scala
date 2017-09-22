package it.unich.jandom.domains.numerical.octagon.fast

import it.unich.jandom.domains.numerical.octagon.variables.VarIndexOps
import it.unich.jandom.domains.numerical.octagon.variables.CountOps._
import it.unich.jandom.domains.numerical.octagon.variables.VarCount
import it.unich.jandom.domains.numerical.octagon.variables.Dimension
import it.unich.jandom.domains.numerical.octagon.variables
import it.unich.jandom.domains.numerical.octagon._
import scalaz._

/**
  * Created by fsestini on 7/11/17.
  *
  * Implementation of matrices that are represented by storing only the lower
  * triangular half of it, as explained in Singh et al.
  */

case class HalfMatrix[A] private[fast] (private[fast] val vec: Vector[A], val dimension: Dimension) {

  def toList: List[A] = grid(dimension).map(p =>  vec(elementIndex(p._1, p._2))).toList
  override def toString: String = {
    val pad = 4
    val maxLength =
      toList
        .map(_.toString.size)
        .max

    (allIndices(dimension)).map(
      (i: Int) =>
      (allIndices(dimension)).map(
        (j: Int) => {
            val el = vec(elementIndex(i,j))
          val res =
            if (i / 2 < j / 2) {
              ""
            } else {
              if (el.toString == "Infinity")
                0x221E.toChar.toString // infty
              else
                el.toString
            }
            (" " * ((maxLength + pad) - res.length)) + res
        }
      ).mkString("")
    ).mkString("\n")
  }

  def this(dimension: Dimension, elem: A) =
    this(variables.Fast.fill(
      variables.Fast.dimToVecSize(dimension)
    )(elem), dimension)

  // The lower part of the DBM, that is, elements at line i, column j,
  // such that i >= j or i = signed(j) [Mine06]
  val lowerIndices: Seq[(Int, Int)] =
    grid(dimension).filter(p => p._1 >= p._2 || p._1 == VarIndexOps.signed(p._2))

  private def elementIndex(i: Int, j: Int): Int =
    if (i < j)
      HalfMatrix.getIndex(j^1, i^1)
    else
      HalfMatrix.getIndex(i, j)

  def update(i: Int, j: Int, x: A): HalfMatrix[A] = {
    require(inDimension(i, j, dimension),
      "HalfMatrix.update(Int,Int,A): index out of bounds")
    new HalfMatrix(vec.updated(elementIndex(i, j), x), dimension)
  }

  def update(updater: (Int, Int) => A): HalfMatrix[A] = {
    val newValues = for (i <- 0 until vec.length)
                       yield (updater.tupled(HalfMatrix.fromIndex(i)))
    new HalfMatrix(newValues.toVector, dimension)
  }

  def apply(i: Int, j: Int): A = {
    require(inDimension(i, j, dimension),
      "Halfmatrix.apply: index out of bounds")
    vec(elementIndex(i, j))
  }

  def combine[B, C](that: HalfMatrix[B], f: (A, B) => C): HalfMatrix[C] = {
    require(this.dimension == that.dimension,
      "HalfMatrix.combine: dimension mismatch")
    val mat = (this.vec zip that.vec) map (f.tupled)
    new HalfMatrix(mat, dimension)
  }

  def toSeq: Seq[A] = vec

}

object HalfMatrix {

  private def getIndex(i: Int, j: Int): Int = j + ((i+1)*(i+1))/2
  private def fromIndex(x: Int): (Int, Int) = {
    def aux(i: Int, j: Int): (Int, Int) = {
      if (getIndex(i, j) == x)
        (i, j)
      else if (getIndex(i + 1, j) <= x)
        aux(i + 1, j)
      else
        aux(i, j + 1)
    }

    aux(0,0)
  }

  def apply[A](f: (Int, Int) => A, nOfVars: VarCount): HalfMatrix[A] = {
    val vec : Vector[A] =
      variables.Fast.allIndices(variables.Fast.dimToVecSize(varCountToDim(nOfVars))).map(
        (k) => {
          val (i,j) = fromIndex(k)
          f(i,j)
        }
      ).toVector
    new HalfMatrix[A] (vec, varCountToDim(nOfVars))
  }
}


object HalfMatrixMatrixInstance {

  implicit val halfMatrixIsMatrix: Matrix[HalfMatrix] = new Matrix[HalfMatrix] {

    def combine[A, B, C](f: (A, B) => C)
      (ma: HalfMatrix[A], mb: HalfMatrix[B]) = ma.combine(mb, f)

    def update[A](i: Int, j: Int, x: A)(m: HalfMatrix[A]) = m.update(i, j, x)

    def update[A](f: (Int, Int) => A)(m: HalfMatrix[A]) = m.update(f)

    def get[A](i: Int, j: Int)(m: HalfMatrix[A]): A = m(i, j)

    def foldMap[A, B](fa: HalfMatrix[A])(f: (A) => B)(implicit F: Monoid[B]): B =
      fa.toSeq.toList.map(f).fold(F.zero)((x, y) => F.append(x, y))

    def foldRight[A, B](fa: HalfMatrix[A], z: => B)(f: (A, => B) => B): B =
      fa.toSeq.toList.foldRight(z)((x, y) => f(x, y))

    def pure[A](dimension: Dimension, x: A) = new HalfMatrix(dimension, x)

    def make[A](f: ((Int, Int) => A), dimension: Dimension) =
      HalfMatrix(f, dimToVarCount(dimension))

    def dimension[A](m: HalfMatrix[A]): Dimension = m.dimension

  }
}
