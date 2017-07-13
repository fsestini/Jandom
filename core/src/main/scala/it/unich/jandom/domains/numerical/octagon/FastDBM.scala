package it.unich.jandom.domains.numerical.octagon

import breeze.numerics.pow
import scala.language.higherKinds

import scalaz.{Applicative, Apply, Monoid}

/**
  * Created by fsestini on 7/10/17.
  */

//case class DenseDBM[M[_], A](m: M[A], nOfVariables: Int, e: Matrix[M])
//  extends FastDBM[M, A] {
//
//  private def signed(i: Int): Int = if (i % 2 != 0) (i + 1) else (i - 1)
//
//  private def strengthen(m: M[A])(implicit ifield: InfField[A]): M[A] = {
//    val updater: (Int, Int) => Option[A] = (i: Int, j: Int) => for {
//      a <- e.get(i, j)(m)
//      b <- e.get(i, signed(i))(m)
//      c <- e.get(signed(j), j)(m)
//    } yield ifield.min(a, ifield.half(ifield.+(b, c)))
//    e.update(updater)(m)
//  }
//
//  def strongClosure() = new DenseDBM(strengthen(close(m)), nOfVariables, e)
//
//  def incrementalClosure() = ???
//
//  def union(that: DenseDBM[M, A])(implicit ifield: InfField[A]) =
//    new DenseDBM[M, A](e.combine(ifield.max)(m, that.m), nOfVariables, e)
//
//  def intersection(that: DenseDBM[M, A])(implicit ifield: InfField[A]): DenseDBM[M, A] =
//    new DenseDBM(e.combine(ifield.min)(m, that.m), nOfVariables, e)
//
//  def all[A](l: List[A], p: A => Boolean): Boolean =
//    l.foldRight(true)((x, y) => p(x) && y)
//
//  def compare(other: DenseDBM[M, A])(implicit pos: Poset[A]): Option[Ordering] = {
//    nOfVariables == other.nOfVariables match {
//      case true => {
//        val qqq: M[Option[Ordering]] = e.combine(pos.compare)(m, other.m)
//        val mm: List[Option[Ordering]] = e.toList(qqq)
//        Applicative[Option].sequence(mm) match {
//          case Some(list) => {
//            (all(list, _ == LT()), all(list, _ == EQ()), all(list, _ == GT())) match {
//              case (true, _, _) => Some(LT())
//              case (_, true, _) => Some(EQ())
//              case (_, _, true) => Some(GT())
//              case _ => None
//            }
//          }
//          case None => None
//        }
//      }
//      case false => None
//    }
//  }
//
//  def update(i: Int, j: Int, x: A): DenseDBM[M, A] =
//    new DenseDBM[M, A](e.update(i, j, x)(m), nOfVariables, e)
//  def update(updater: (Int, Int) => Option[A]): DenseDBM[M, A] =
//    new DenseDBM[M, A](e.update(updater)(m), nOfVariables, e)
//  def get(i: Int, j: Int): Option[A] = e.get(i,j)(m)
//  def combine[B, C](f: (A, B) => C, other: DenseDBM[M, B]): DenseDBM[M, C] =
//    new DenseDBM(e.combine(f)(m, other.m), nOfVariables, e)
//}

// Sparsity D = 1 - (nni / (2n^2 + 2n))
// Switching between DBMs: the sparsity can increase, for instance during
// widening. Recovering sparsity information and  independent components has a
// quadratic worst case complexity, so we only perform it by piggybacking on
// the closure operator. We also use closure computations as switching points.

sealed trait FastDBM[M[_], A]
case class DenseDBM[M[_], A](m: M[A], rdbm: RawDBM[M]) extends FastDBM[M, A] { }
case class SparseDBM[M[_], A](m: M[A], rdbm: RawDBM[M]) extends FastDBM[M, A] { }

// We store the independent components as a linked list of linked lists of
// variable indices.
// Independent components correspond to submatrices in the complete DBM,
// and these can be dense or sparse. The exact sparsity is computed
// on-the-fly before closure.

// Octagon operators on the decomposed type:
// Apart from closure (which is specialized according to sparsity), we apply the
// standard operators independently on each submatrix.

// NNI for decomposed DBMs is computed as the sum of its submatrices' NNI.
case class DecomposedDBM[M[_], A](completeDBM: M[A],
                                  indepComponents: List[List[Int]],
                                  evidence: RawDBM[M])
  extends FastDBM[M, A] { }

// The Top type can be seen as a degenerate case of the decomposed one with
// an empty set of independent components. Thus, we can reuse the operators of
// the decomposed type.
case class TopDBM[M[_], A]() extends FastDBM[M, A] { }
case class BottomDBM[M[_], A]() extends FastDBM[M, A] { }

object FastDBMTypeclasses {

  def fastDBMIsDBM[M[_]](implicit e: Matrix[M])
  : DifferenceBoundMatrix[({type T[A] = FastDBM[M, A] })#T] = ???
//    new DifferenceBoundMatrix[({type T[A] = FastDBM[M, A] })#T] {
//
//      def strongClosure[A](m: FastDBM[M, A])(implicit evidence: InfField[A]): FastDBM[M, A] = ???
//      def incrementalClosure[A](m: FastDBM[M, A])(implicit evidence: InfField[A]): FastDBM[M, A] = ???
//
//      def bottomDBM[A]: FastDBM[M, A] = BottomDBM()
//
//      def get[A](i: Int, j: Int)(m: FastDBM[M, A]): Option[A] = ???
//
//      // join and meet must act on submatrices if the DBM is decomposed!
//      def union[A](x: FastDBM[M, A], y: FastDBM[M, A]): FastDBM[M, A] = ???
//      def intersection[A](x: FastDBM[M, A], y: FastDBM[M, A]): FastDBM[M, A] = ???
//    }

}

//object DenseDBMIsDBM {
//  def denseDBMIsDBM[M[_]](implicit me: Matrix[M]): DifferenceBoundMatrix[({type L[a] = DenseDBM[M, a]})#L] =
//    new DifferenceBoundMatrix[({type L[a] = DenseDBM[M, a]})#L] {
//      def compare[A](x: DenseDBM[M, A], y: DenseDBM[M, A])(implicit evidence: PosetConstraint[A]): Option[Ordering] =
//        x.compare(y)
//      def update[A](i: Int, j: Int, x: A)(dbm: DenseDBM[M, A]): DenseDBM[M, A] =
//        dbm.update(i, j, x)
//      def update[A](updater: (Int, Int) => Option[A])(m: DenseDBM[M, A]): DenseDBM[M, A] =
//        m.update(updater)
//
//      def get[A](i: Int, j: Int)(m: DenseDBM[M, A]): Option[A] = m.get(i, j)
//      def combine[A, B, C](f: (A, B) => C)(ma: DenseDBM[M, A], mb: DenseDBM[M, B]): DenseDBM[M, C] =
//        ma.combine(f, mb)
//
//      def foldMap[A, B](fa: DenseDBM[M, A])(f: (A) => B)(implicit F: Monoid[B]): B = fa match {
//        case DenseDBM(m, _, _) => me.foldMap(m)(f)
//      }
//      def foldRight[A, B](fa: DenseDBM[M, A], z: => B)(f: (A, B) => B): B = fa match {
//        case DenseDBM(m, _, _) => me.foldRight(m, z)(f)
//      }
//      def map[A, B](fa: DenseDBM[M, A])(f: (A) => B): DenseDBM[M, B] = fa match {
//        case DenseDBM(m, n, e) => DenseDBM(me.map(m)(f), n, e)
//      }
//
//    }
//}

//
//object DoubleFunDense {
//  def apply(nOfVariables: Int)(implicit evidence: Matrix[FunMatrix]): DenseDBM[FunMatrix, Double] =
//    new DenseDBM[FunMatrix, Double](
//      new FunMatrix[Double](
//        (_,_) => Some(Double.PositiveInfinity),
//        2 * nOfVariables),
//      nOfVariables, evidence)
//}
