package it.unich.jandom.domains.numerical.octagon

import scalaz.{Applicative, Apply, Monoid}
import scalaz.std.option._, scalaz.std.list._

/**
  * Created by fsestini on 7/13/17.
  */

// FunMatrix-based raw DBM implementation
case class FunRawDBM[A](nOfVars: Int, m: FunMatrix[A]) { }

object FunRawDBM {
  def apply[A](nOfVars: Int, fun: (Int, Int) => A) =
    new FunRawDBM(nOfVars, new FunMatrix(fun, nOfVars * 2))
}

object FunDBMIsRawDBM {

  def funDBMIsRawDBM(implicit e: Matrix[FunMatrix]): RawDBM[FunRawDBM] = new RawDBM[FunRawDBM] {

    private def cloTransform[A](k: Int)(m: FunRawDBM[A])(implicit ifield: InfField[A]): FunRawDBM[A] = {

      val f: (Int, Int) => A = (i: Int, j: Int) => {
        val l: List[Option[A]] = List(
          get(i,j)(m),
          Apply[Option].apply2(get(i, 2 * k - 1)(m), get(2 * k - 1, j)(m))(ifield.+),
          Apply[Option].apply2(get(i, 2 * k)(m), get(2 * k, j)(m))(ifield.+),
          Apply[Option].apply2(
            get(i, 2 * k - 1)(m),
            Apply[Option].apply2(
              get(2 * k - 1, 2 * k)(m),
              get(2 * k, j)(m))(ifield.+)
          )(ifield.+),
          Apply[Option].apply2(
            get(i, 2 * k)(m),
            Apply[Option].apply2(
              get(2 * k, 2 * k - 1)(m),
              get(2 * k - 1, j)(m))(ifield.+)
          )(ifield.+)
        )
        Applicative[Option].sequence(l) match {
          case Some(resList) => resList.fold(ifield.infinity)(ifield.min)
          case None => throw new IllegalArgumentException()
        }
      }
      update(f)(m)
    }

    private def close[A](m: FunRawDBM[A])(implicit ifield: InfField[A]): FunRawDBM[A] = {
      var x = m
      for (k <- 1 to m.nOfVars) {
        x = cloTransform(k)(x)
      }
      x
    }

    private def strengthen[A](m: FunRawDBM[A])(implicit ifield: InfField[A]): FunRawDBM[A] = {
      val updater: (Int, Int) => A = (i: Int, j: Int) => {
        val o = for {
          a <- get(i, j)(m)
          b <- get(i, signed(i))(m)
          c <- get(signed(j), j)(m)
        } yield ifield.min(a, ifield.half(ifield.+(b, c)))
        o match {
          case Some(r) => r
          case None => throw new IllegalArgumentException()
        }
      }
      update(updater)(m)
    }

    private def signed(i: Int): Int = if (i % 2 != 0) i + 1 else i - 1

    def nOfVars[A](m: FunRawDBM[A]): VarIndex =
      m match { case FunRawDBM(n, _) => n }

    def denseStrongClosure[A](m: FunRawDBM[A])(implicit e: InfField[A])
    : (FunRawDBM[A], VarIndex, List[List[VarIndex]]) =
      (strengthen(close(m)), ???, ???)

    def denseStrongClosure[A](m: FunRawDBM[A], indices: Seq[VarIndex])(implicit e: InfField[A])
    : (FunRawDBM[A], VarIndex, List[List[VarIndex]]) = ???

    def sparseStrongClosure[A](m: FunRawDBM[A])(implicit e: InfField[A])
    : (FunRawDBM[A], VarIndex, List[List[VarIndex]]) = ???

    def sparseStrongClosure[A](m: FunRawDBM[A], indices: Seq[VarIndex])(implicit e: InfField[A])
    : (FunRawDBM[A], VarIndex, List[List[VarIndex]]) = ???

    def denseIncrementalClosure[A](m: FunRawDBM[A])(implicit e: InfField[A])
    : (FunRawDBM[A], VarIndex, List[List[VarIndex]]) = ???

    def sparseIncrementalClosure[A](m: FunRawDBM[A])(implicit e: InfField[A])
    : (FunRawDBM[A], VarIndex, List[List[VarIndex]]) = ???

    def combine[A, B, C](f: (A, B) => C, xScope: Seq[VarIndex], yScope: Seq[VarIndex])(ma: FunRawDBM[A], mb: FunRawDBM[B]): FunRawDBM[C] = ???

    def combine[A, B, C](f: (A, B) => C)(ma: FunRawDBM[A], mb: FunRawDBM[B]): FunRawDBM[C] =
      (ma, mb) match {
        case (FunRawDBM(n1, m1), FunRawDBM(n2, m2)) => FunRawDBM(n1, e.combine(f)(m1, m2))
      }

    def update[A](i: Int, j: Int, x: A)(m: FunRawDBM[A]): FunRawDBM[A] =
      m match { case FunRawDBM(n, mx) => FunRawDBM(n, e.update(i, j, x)(mx)) }

    def update[A](updater: (Int, Int) => A)(m: FunRawDBM[A]): FunRawDBM[A] =
      m match { case FunRawDBM(n, mx) => FunRawDBM(n, e.update(updater)(mx)) }

    def get[A](i: Int, j: Int)(m: FunRawDBM[A]): Option[A] =
      m match { case FunRawDBM(n, mx) => e.get(i, j)(mx) }

    def foldMap[A, B](fa: FunRawDBM[A])(f: (A) => B)(implicit F: Monoid[B]): B =
      fa match { case FunRawDBM(n, mx) => e.foldMap(mx)(f) }

    def foldRight[A, B](fa: FunRawDBM[A], z: => B)(f: (A, => B) => B): B =
      fa match { case FunRawDBM(n, mx) => e.foldRight(mx, z)(f) }

  }
}