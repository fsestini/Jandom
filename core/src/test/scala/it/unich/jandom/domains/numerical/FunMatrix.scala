package it.unich.jandom.domains.numerical.octagon
import scalaz._
import it.unich.jandom.domains.numerical.octagon.variables.CountOps
import it.unich.jandom.domains.numerical.octagon.variables.VarCount
import it.unich.jandom.domains.numerical.octagon.variables.Dimension


// Simple implementation of a functional square matrix.
// For test purposes only.
class FunMatrix[A](private val fun: (Int, Int) => A, val dimension: Dimension) {

  def ==(that: FunMatrix[A]): Boolean = {
    if (dimension != that.dimension)
      false
    else {
      val thisTupled = this.fun.tupled
      val thatTupled = that.fun.tupled
      CountOps.grid(dimension).forall(idx => thisTupled(idx) == thatTupled(idx))
    }
  }

  def update(i: Int, j: Int, x: A): FunMatrix[A] = {
    require(CountOps.inDimension(i, j, dimension),
      "Can't update (" + i + "," + j + "), dimension is " + dimension)
    update((ii, jj) => if (ii == i && jj == j) x else fun(ii, jj))
  }

  def update(updater: (Int, Int) => A): FunMatrix[A] =
    new FunMatrix(updater, dimension)

  def apply(i: Int, j: Int): A = {
    require(CountOps.inDimension(i, j, dimension),
      "Can't apply (" + i + "," + j + "), dimension is " + dimension)
    fun(i, j)
  }

  def combine[B, C](that: FunMatrix[B], f: (A, B) => C): FunMatrix[C] = {
    val newFun: (Int, Int) => C = (i, j) => f(this(i, j), that(i, j))
    new FunMatrix[C](newFun, dimension)
  }

  private def filterOut(l: List[Option[A]]): List[A] = l match {
    case (x :: xs) => x match {
      case Some(y) => y :: filterOut(xs)
      case None => filterOut(xs)
    }
    case _ => Nil
  }

  def toList: List[A] = CountOps.grid(dimension).map(p => fun(p._1, p._2)).toList

  override def toString: String = {
    val pad = 4
    val maxLength =
      toList
        .map(_.toString.size)
        .max

    (CountOps.allIndices(dimension)).map(
      (i: Int) =>
      (CountOps.allIndices(dimension)).map(
        (j: Int) => {
          val res =
            if (fun(i,j) == Double.PositiveInfinity)
              0x221E.toChar.toString // infty
            else
              fun(i,j).toString
          (" " * ((maxLength + pad) - res.length)) + res
        }
      ).mkString("")
    ).mkString("\n")
  }
}

object FunMatrix {
  def apply[A](fun: (Int, Int) => A, dimension: Dimension) =
    new FunMatrix(fun, dimension)
}

object FunMatrixMatrixInstance {
  implicit val funMatrixIsMatrix: Matrix[FunMatrix] = new Matrix[FunMatrix] {
    def combine[A, B, C](f: (A, B) => C)
                        (ma: FunMatrix[A], mb: FunMatrix[B]): FunMatrix[C] =
      ma.combine(mb, f)
    def update[A](i: Int, j: Int, x: A)(m: FunMatrix[A]): FunMatrix[A] =
      m.update(i, j, x)
    def update[A](updater: (Int, Int) => A)(m: FunMatrix[A]): FunMatrix[A] =
      m.update(updater)
    def get[A](i: Int, j: Int)(m: FunMatrix[A]): A = m(i, j)
    def foldMap[A, B](fa: FunMatrix[A])(f: (A) => B)(implicit F: Monoid[B]): B =
      fa.toList.map(f).fold(F.zero)((x, y) => F.append(x, y))

    def foldRight[A, B](fa: FunMatrix[A], z: => B)(f: (A, => B) => B): B =
      fa.toList.foldRight(z)((x, y) => f(x, y))

    def pure[A](dimension: Dimension, x: A): FunMatrix[A] =
      new FunMatrix[A]((_, _) => x, dimension)
    def make[A](f: (Int, Int) => A, dimension: Dimension): FunMatrix[A] = FunMatrix(f, dimension)
    def dimension[A](m: FunMatrix[A]): Dimension = m.dimension
  }
}
