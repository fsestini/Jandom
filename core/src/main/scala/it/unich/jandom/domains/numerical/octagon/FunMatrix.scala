package it.unich.jandom.domains.numerical.octagon
import scalaz.Monoid

case class FunMatrix[A](fun: (Int, Int) => Option[A], dimension: Int) {
  def update(i: Int, j: Int, x: A): FunMatrix[A] =
    update((ii: Int, jj: Int) => (ii == i, jj == j) match {
      case (true, true) => Some(x)
      case _ => fun(i,j)
    })

  def update(updater: (Int, Int) => Option[A]): FunMatrix[A] =
    new FunMatrix(updater, dimension)

  def apply(i: Int, j: Int): Option[A] =
    (1 <= i, i <= dimension, 1 <= j, j <= dimension) match {
      case (true, true, true, true) => fun(i, j)
      case _ => None
    }

  def combine[B, C](other: FunMatrix[B], f: (A, B) => C): FunMatrix[C] = {
    val newFun: (Int, Int) => Option[C] = (i: Int, j: Int) => for {
      x <- this(i, j)
      y <- other(i, j)
    } yield f(x, y)
    new FunMatrix[C](newFun, dimension)
  }

  private def filterOut[A](l: List[Option[A]]): List[A] = l match {
    case (x :: xs) => x match {
      case Some(y) => y :: filterOut(xs)
      case None => filterOut(xs)
    }
    case _ => Nil
  }

  def toList: List[A] = {
    val indexes: List[(Int, Int)] = (for {
      x <- 1 to dimension
      y <- 1 to dimension}
      yield (x, y)).toList
    val l: List[Option[A]] = for ((i, j) <- indexes) yield fun(i, j)
    filterOut(l)
  }
}

object FunMatrixIsMatrix {
  implicit val FunMatrixIsMatrix: Matrix[FunMatrix] = new Matrix[FunMatrix] {
    def update[A](i: Int, j: Int, x: A)(m: FunMatrix[A]): FunMatrix[A] =
      m.update(i, j, x)
    def update[A](updater: (Int, Int) => Option[A])(m: FunMatrix[A]): FunMatrix[A] =
      m.update(updater)
    def get[A](i: Int, j: Int)(m: FunMatrix[A]): Option[A] = m(i,j)
    def combine[A, B, C](f: (A, B) => C)(ma: FunMatrix[A], mb: FunMatrix[B]): FunMatrix[C] =
      ma.combine(mb, f)
    // def empty[A](dim: Int): FunMatrix[A] = new FunMatrix[A]((_, _) => None, dim)

    def foldMap[A, B](fa: FunMatrix[A])(f: (A) => B)(implicit F: Monoid[B]): B =
      this.foldRight(fa, F.zero)((x, y) => F.append(f(x), y))

    def foldRight[A, B](fa: FunMatrix[A], z: => B)(f: (A, B) => B): B =
      fa.toList.foldRight(z)(f)

    def map[A, B](fa: FunMatrix[A])(f: (A) => B): FunMatrix[B] = {
      val updater: (Int, Int) => Option[B] = (i: Int, j: Int) => for {
        x <- fa(i, j)
      } yield f(x)
      new FunMatrix[B](updater, fa.dimension)
    }

  }
}

object FunMatrix {
  def apply[A](dim: Int, init: A) = {
    new FunMatrix((_,_) => Some(init), dim)
  }
}
