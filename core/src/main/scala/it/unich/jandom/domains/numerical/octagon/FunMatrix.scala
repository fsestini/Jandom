package it.unich.jandom.domains.numerical.octagon
import scalaz._

import scala.language.higherKinds

// Simple implementation of a functional square matrix.
// For test purposes only.
case class FunMatrix[A](fun: (Int, Int) => A, dimension: Int) {

  def update(i: Int, j: Int, x: A): FunMatrix[A] =
    update((ii,jj) => if ((ii,jj) == (i, j)) x else fun(i,j))

  def update(updater: (Int, Int) => A): FunMatrix[A] =
    new FunMatrix(updater, dimension)

  def apply(i: Int, j: Int): Option[A] =
    (1 <= i, i <= dimension, 1 <= j, j <= dimension) match {
      case (true, true, true, true) => Some(fun(i, j))
      case _ => None
    }

  def combine[B, C](other: FunMatrix[B], f: (A, B) => C): FunMatrix[C] = {
    val newFun: (Int, Int) => C = (i: Int, j: Int) => {
      val res = for {
        x <- this(i, j)
        y <- other(i, j)
      } yield f(x, y)
      res match {
        case Some(r) => r
        case None => throw new IllegalArgumentException()
      }
    }
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
      x <- 1 to dimension
      y <- 1 to dimension
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
    def get[A](i: Int, j: Int)(m: FunMatrix[A]): Option[A] = m(i, j)
    def foldMap[A, B](fa: FunMatrix[A])(f: (A) => B)(implicit F: Monoid[B]): B =
      fa.toList.map(f).fold(F.zero)((x, y) => F.append(x, y))

    def foldRight[A, B](fa: FunMatrix[A], z: => B)(f: (A, => B) => B): B =
      fa.toList.foldRight(z)((x, y) => f(x, y))
  }
}

