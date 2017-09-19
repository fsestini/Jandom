package it.unich.jandom.domains.numerical.octagon.fast
import it.unich.jandom.domains.numerical.octagon.variables.VarIndexOps
import it.unich.jandom.domains.numerical.octagon.variables.CountOps._
import it.unich.jandom.domains.numerical.octagon.variables.VarCount
import it.unich.jandom.domains.numerical.octagon.variables.Dimension
import it.unich.jandom.domains.numerical.octagon.variables

/**
  * Created by fsestini on 7/11/17.
  *
  * Implementation of matrices that are represented by storing only the lower
  * triangular half of it, as explained in Singh et al.
  */

class HalfMatrix[A] private[fast] (private[fast] val vec: Vector[A], val dimension: Dimension) {

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

  private def elementIndex(i: Int, j: Int): Int =
    if (i < j)
      getIndex(j^1, i^1)
    else
      getIndex(i, j)

  def update(i: Int, j: Int, x: A): HalfMatrix[A] = {
    require(inDimension(i, j, dimension))
    new HalfMatrix(vec.updated(elementIndex(i, j), x), dimension)
  }

  def update(updater: (Int, Int) => A): HalfMatrix[A] = {
    val newValues = for (i <- 0 until vec.length)
                       yield (updater.tupled(fromIndex(i)))
    new HalfMatrix(newValues.toVector, dimension)
  }

  def apply(i: Int, j: Int): A = {
    require(inDimension(i, j, dimension))
    vec(elementIndex(i, j))
  }

  def combine[B, C](that: HalfMatrix[B], f: (A, B) => C): HalfMatrix[C] = {
    require(this.dimension == that.dimension)
    val mat = (this.vec zip that.vec) map (f.tupled)
    new HalfMatrix(mat, dimension)
  }

  def toSeq: Seq[A] = vec

}

object HalfMatrix {
  def apply[A](f: (Int, Int) => A, nOfVars: VarCount): HalfMatrix[A] = {
    val vec : Vector[A] = (allIndices(varCountToDim(nOfVars))).map(
      (i) => (allIndices(varCountToDim(nOfVars))).map(
        (j) =>
        f(i,j)
      )
    ).flatten.toVector
    new HalfMatrix[A] (vec, varCountToDim(nOfVars))
  }
}
