package it.unich.jandom.domains.numerical.octagon
import scalaz.Monoid

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

// FunMatrix-based raw DBM implementation
case class FunDBM[A](nOfVars: Int, m: FunMatrix[A]) { }

object FunDBM {
  def apply(nOfVars: Int, fun: (Int, Int) => A) =
    new FunDBM(nOfVars, new FunMatrix(fun, nOfVars * 2))
}

object FunDBMIsRawDBM {
  implicit val funDBMIsRawDBM: RawDBM[FunDBM] = new RawDBM[FunDBM] {

    private def cloTransform[M[_]](k: Int)(m: M[Double])
      (implicit e: Matrix[M], ifield: InfField[A]): M[Double] = {

      val f: (Int, Int) => A = (i: Int, j: Int) => {
        val l: List[Option[A]] = List(
          e.get(i,j)(m),
          Apply[Option].apply2(e.get(i, 2 * k - 1)(m), e.get(2 * k - 1, j)(m))(_ + _),
          Apply[Option].apply2(e.get(i, 2 * k)(m), e.get(2 * k, j)(m))(_ + _),
          Apply[Option].apply2(
            e.get(i, 2 * k - 1)(m),
            Apply[Option].apply2(
              e.get(2 * k - 1, 2 * k)(m),
              e.get(2 * k, j)(m))(_ + _)
          )(_ + _),
          Apply[Option].apply2(
            e.get(i, 2 * k)(m),
            Apply[Option].apply2(
              e.get(2 * k, 2 * k - 1)(m),
              e.get(2 * k - 1, j)(m))(_ + _)
          )(_ + _)
        )
        Applicative[Option].sequence(l) match {
          case Some(resList) => Some(resList.fold(ifield.infinity)(ifield.min _))
          case None => None
        }
      }
      e.update(f)(m)
    }

    private def close[M[_]](m: M[Double])(implicit me: Matrix[M]): M[Double] = {
      var x: M[Double] = m
      for (k <- 1 to nOfVariables) {
        x = cloTransform(k)(x)
      }
      x
    }

    private def strengthen[M[_]](m: M[A])
      (implicit e: InfField[A], me: Matrix[M]): M[A] = {
      val updater: (Int, Int) => Double = (i: Int, j: Int) => for {
        a <- e.get(i, j)(m)
        b <- e.get(i, signed(i))(m)
        c <- e.get(signed(j), j)(m)
      } yield math.min(a, (b + c) / 2)
      e.update(updater)(m)
    }

    private def signed(i: Int): Int = if (i % 2 != 0) (i + 1) else (i - 1)

    def denseStrongClosure[A](m: FunDBM[A])(implicit e: InfField[A]): ClosureRes[A] =
      (strengthen(close(m)), ???, ???)
    def combine[A, B, C](f: (A, B) => C)(ma: FunDBM[A], mb: FunDBM[B]): FunDBM[C] =
      ma.combine(bm, f)
    def update[A](i: VarIndex, j: VarIndex, x: A)(m: FunDBM[A]): FunDBM[A] =
      m.update(i, j, x)
    def update[A](updater: (VarIndex, VarIndex) => Option[A])(m: FunDBM[A]): FunDBM[A] =
      m.update(updater)
    def get[A](i: VarIndex, j: VarIndex)(m: FunDBM[A]): Option[A] = m(i, j)

    def denseStrongClosure[A](m: FunDBM[A], indices: Seq[VarIndex])(implicit e: InfField[A]): ClosureRes[A] = ???
    def combine[A, B, C](f: (A, B) => C, xScope: Seq[VarIndex], yScope: Seq[VarIndex])(ma: FunDBM[A], mb: FunDBM[B]): FunDBM[C] = ???

    def foldMap[A, B](fa: FunDBM[A])(f: (A) => B)(implicit F: Monoid[B]): B = ???
    def foldRight[A, B](fa: FunDBM[A], z: => B)(f: (A, B) => B): B = ???

    def sparseStrongClosure[A](m: FunDBM[A])(implicit e: InfField[A]): ClosureRes[A] = ???
    def sparseStrongClosure[A](m: FunDBM[A], indices: Seq[VarIndex])(implicit e: InfField[A]): ClosureRes[A] = ???
    def denseIncrementalClosure[A](m: FunDBM[A])(implicit e: InfField[A]): ClosureRes[A] = ???
    def sparseIncrementalClosure[A](m: FunDBM[A])(implicit e: InfField[A]): ClosureRes[A] = ???
  }
}

// object FunMatrixIsMatrix {
//   implicit val FunMatrixIsMatrix: Matrix[FunMatrix] = new Matrix[FunMatrix] {
//     def update[A](i: Int, j: Int, x: A)(m: FunMatrix[A]): FunMatrix[A] =
//       m.update(i, j, x)
//     def update[A](updater: (Int, Int) => Option[A])(m: FunMatrix[A]): FunMatrix[A] =
//       m.update(updater)
//     def get[A](i: Int, j: Int)(m: FunMatrix[A]): Option[A] = m(i,j)
//     def combine[A, B, C](f: (A, B) => C)(ma: FunMatrix[A], mb: FunMatrix[B]): FunMatrix[C] =
//       ma.combine(mb, f)
//     // def empty[A](dim: Int): FunMatrix[A] = new FunMatrix[A]((_, _) => None, dim)

//     def foldMap[A, B](fa: FunMatrix[A])(f: (A) => B)(implicit F: Monoid[B]): B =
//       this.foldRight(fa, F.zero)((x, y) => F.append(f(x), y))

//     def foldRight[A, B](fa: FunMatrix[A], z: => B)(f: (A, B) => B): B =
//       fa.toList.foldRight(z)(f)

//     def map[A, B](fa: FunMatrix[A])(f: (A) => B): FunMatrix[B] = {
//       val updater: (Int, Int) => Option[B] = (i: Int, j: Int) => for {
//         x <- fa(i, j)
//       } yield f(x)
//       new FunMatrix[B](updater, fa.dimension)
//     }

//   }
// }

// case class FunMatrix[A](fun: (Int, Int) => Option[A], dimension: Int) {

//   def update(i: Int, j: Int, x: A): FunMatrix[A] =
//     update((ii,jj) => if ((ii,jj) == (i, j)) Some(x) else fun(i,j))

//   def update(updater: (Int, Int) => Option[A]): FunMatrix[A] =
//     new FunMatrix(updater, dimension)

//   def apply(i: Int, j: Int): Option[A] =
//     (1 <= i, i <= dimension, 1 <= j, j <= dimension) match {
//       case (true, true, true, true) => fun(i, j)
//       case _ => None
//     }

//   def combine[B, C](other: FunMatrix[B], f: (A, B) => C): FunMatrix[C] = {
//     val newFun: (Int, Int) => Option[C] = (i: Int, j: Int) => for {
//       x <- this(i, j)
//       y <- other(i, j)
//     } yield f(x, y)
//     new FunMatrix[C](newFun, dimension)
//   }

//   private def filterOut[A](l: List[Option[A]]): List[A] = l match {
//     case (x :: xs) => x match {
//       case Some(y) => y :: filterOut(xs)
//       case None => filterOut(xs)
//     }
//     case _ => Nil
//   }

//   def toList: List[A] = {
//     val indexes: List[(Int, Int)] = (for {
//       x <- 1 to dimension
//       y <- 1 to dimension}
//       yield (x, y)).toList
//     val l: List[Option[A]] = for ((i, j) <- indexes) yield fun(i, j)
//     filterOut(l)
//   }
// }

// object FunMatrixIsMatrix {
//   implicit val FunMatrixIsMatrix: Matrix[FunMatrix] = new Matrix[FunMatrix] {
//     def update[A](i: Int, j: Int, x: A)(m: FunMatrix[A]): FunMatrix[A] =
//       m.update(i, j, x)
//     def update[A](updater: (Int, Int) => Option[A])(m: FunMatrix[A]): FunMatrix[A] =
//       m.update(updater)
//     def get[A](i: Int, j: Int)(m: FunMatrix[A]): Option[A] = m(i,j)
//     def combine[A, B, C](f: (A, B) => C)(ma: FunMatrix[A], mb: FunMatrix[B]): FunMatrix[C] =
//       ma.combine(mb, f)
//     // def empty[A](dim: Int): FunMatrix[A] = new FunMatrix[A]((_, _) => None, dim)

//     def foldMap[A, B](fa: FunMatrix[A])(f: (A) => B)(implicit F: Monoid[B]): B =
//       this.foldRight(fa, F.zero)((x, y) => F.append(f(x), y))

//     def foldRight[A, B](fa: FunMatrix[A], z: => B)(f: (A, B) => B): B =
//       fa.toList.foldRight(z)(f)

//     def map[A, B](fa: FunMatrix[A])(f: (A) => B): FunMatrix[B] = {
//       val updater: (Int, Int) => Option[B] = (i: Int, j: Int) => for {
//         x <- fa(i, j)
//       } yield f(x)
//       new FunMatrix[B](updater, fa.dimension)
//     }

//   }
// }

// object FunMatrix {
//   def apply[A](dim: Int, init: A) = {
//     new FunMatrix((_,_) => Some(init), dim)
//   }
// }
