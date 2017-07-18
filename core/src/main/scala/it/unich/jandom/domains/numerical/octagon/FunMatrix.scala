package it.unich.jandom.domains.numerical.octagon
import scalaz._

import scala.language.higherKinds

// Simple implementation of a functional square matrix.
// For test purposes only.
case class FunMatrix[A](fun: (Int, Int) => A, dimension: Int) {

  def update(i: Int, j: Int, x: A): FunMatrix[A] = {
    require(0 <= i && i < dimension && 0 <= j && j < dimension)
    update((ii, jj) => if ((ii, jj) == (i, j)) x else fun(i, j))
  }

  def update(updater: (Int, Int) => A): FunMatrix[A] =
    new FunMatrix(updater, dimension)

  def apply(i: Int, j: Int): A = {
    require(0 <= i && i < dimension && 0 <= j && j < dimension)
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

  def toList: List[A] = {
    val indexes: List[(Int, Int)] = (for {
      x <- 0 until dimension
      y <- 0 until dimension
    } yield (x, y)).toList
    for ((i, j) <- indexes) yield fun(i, j)
  }
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

    def pure[A](dimension: Int, x: A): FunMatrix[A] =
      new FunMatrix[A]((_, _) => x, dimension)
  }
}

