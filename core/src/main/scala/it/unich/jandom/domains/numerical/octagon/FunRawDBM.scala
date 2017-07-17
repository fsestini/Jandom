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

  def funDBMIsRawDBM(implicit e: Matrix[FunMatrix]): DenseSparseDBM[FunRawDBM] =
    new DenseSparseDBM[FunRawDBM] {

    private def cloTransform[A](k: Int)(m: FunRawDBM[A])(implicit ifield: InfField[A]): FunRawDBM[A] = {

      val f: (Int, Int) => A = (i: Int, j: Int) => {
        val l: List[A] = List(
          get(i,j)(m),
          ifield.+(
            get(i, 2 * k - 1)(m),
            get(2 * k - 1, j)(m)
          ),
          ifield.+(
            get(i, 2 * k)(m),
            get(2 * k, j)(m)
          ),
          ifield.+(
            get(i, 2 * k - 1)(m),
            ifield.+(
              get(2 * k - 1, 2 * k)(m),
              get(2 * k, j)(m))
          ),
          ifield.+(
            get(i, 2 * k)(m),
            ifield.+(
              get(2 * k, 2 * k - 1)(m),
              get(2 * k - 1, j)(m))
          )
        )

        l.fold(ifield.infinity)(ifield.min)
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
        val a = get(i, j)(m)
        val b = get(i, signed(i))(m)
        val c = get(signed(j), j)(m)
        ifield.min(a, ifield.half(ifield.+(b, c)))
      }
      update(updater)(m)
    }

    private def signed(i: Int): Int = if (i % 2 != 0) i + 1 else i - 1

    def update[A](i: Int, j: Int, x: A)(m: FunRawDBM[A]): FunRawDBM[A] =
      m match { case FunRawDBM(n, mx) => FunRawDBM(n, e.update(i, j, x)(mx)) }

    def update[A](updater: (Int, Int) => A)(m: FunRawDBM[A]): FunRawDBM[A] =
      m match { case FunRawDBM(n, mx) => FunRawDBM(n, e.update(updater)(mx)) }

    def get[A](i: Int, j: Int)(m: FunRawDBM[A]): A =
      m match { case FunRawDBM(n, mx) => e.get(i, j)(mx) }

    def foldMap[A, B](fa: FunRawDBM[A])(f: (A) => B)(implicit F: Monoid[B]): B =
      fa match { case FunRawDBM(n, mx) => e.foldMap(mx)(f) }

    def foldRight[A, B](fa: FunRawDBM[A], z: => B)(f: (A, => B) => B): B =
      fa match { case FunRawDBM(n, mx) => e.foldRight(mx, z)(f) }

    def extract[A](is: Seq[it.unich.jandom.domains.numerical.octagon.VarIndex])(m: it.unich.jandom.domains.numerical.octagon.FunRawDBM[A]): it.unich.jandom.domains.numerical.octagon.FunRawDBM[A] = ???
    def forget[A](v: it.unich.jandom.domains.numerical.octagon.VarIndex)(m: it.unich.jandom.domains.numerical.octagon.FunRawDBM[A]): it.unich.jandom.domains.numerical.octagon.FunRawDBM[A] = ???
    def incrementalClosure[A](v: it.unich.jandom.domains.numerical.octagon.VarIndex)(m: it.unich.jandom.domains.numerical.octagon.FunRawDBM[A])(implicit e: it.unich.jandom.domains.numerical.octagon.InfField[A]): Option[(it.unich.jandom.domains.numerical.octagon.FunRawDBM[A], it.unich.jandom.domains.numerical.octagon.NNI, List[List[it.unich.jandom.domains.numerical.octagon.VarIndex]])] = ???
    def pour[A](source: it.unich.jandom.domains.numerical.octagon.FunRawDBM[A])(dest: it.unich.jandom.domains.numerical.octagon.FunRawDBM[A]): it.unich.jandom.domains.numerical.octagon.FunRawDBM[A] = ???
    def strongClosure[A](m: it.unich.jandom.domains.numerical.octagon.FunRawDBM[A])(implicit e: it.unich.jandom.domains.numerical.octagon.InfField[A]): Option[(it.unich.jandom.domains.numerical.octagon.FunRawDBM[A], it.unich.jandom.domains.numerical.octagon.NNI, List[List[it.unich.jandom.domains.numerical.octagon.VarIndex]])] = ???
    def varIndices[A](m: it.unich.jandom.domains.numerical.octagon.FunRawDBM[A]): Seq[it.unich.jandom.domains.numerical.octagon.VarIndex] = ???

  }
}
