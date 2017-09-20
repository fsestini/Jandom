package it.unich.jandom.domains.numerical.octagon.fast

import breeze.numerics.pow
import it.unich.jandom.domains.numerical.octagon._
import variables._
import CountOps._
import scala.language.higherKinds
import scalaz.{Applicative, Apply, Monoid, Traverse}
import scalaz.std.option._
import scalaz.std.list._

// Sparsity D = 1 - (nni / (2n^2 + 2n))
// Switching between DBMs: the sparsity can increase, for instance during
// widening. Recovering sparsity information and  independent components has a
// quadratic worst case complexity, so we only perform it by piggybacking on
// the closure operator. We also use closure computations as switching points.

// DifferenceBoundMatrix instance for CFastDBM.
object CFDBMInstance {

  def instance[M[_], SM[_]](implicit mev: MEvidence[M, SM]) =
    new DifferenceBoundMatrix[({ type T[S, A] = CFastDBM[M, SM, S, A] })#T] {

      type PosetConstraint[A] = InfField[A]

      def update[S <: DBMState, A](f: (Int, Int) => A)(dbm: CFastDBM[M, SM, S, A])
        (implicit ifield: InfField[A]): ExistsM[A] =
        // Applying an arbitrary map destroys every knowledge about relations
        // between variables, so we just return DBM with unknown closure state
        // (because it could be either non-closed or bottom, which is closed).
        cfastInnerMatrix(dbm).map(mev.ds.update(f)) match {
          case Some(m) => packEx(NCFast(FullDBM(m, mev)))
          case None => packEx(BottomFast(FastDBMUtils.nOfVars(dbm)))
        }

      def incrementalClosure[S <: DBMState, A](v: VarIndex)
                               (dbm: CFastDBM[M, SM, S, A])
        (implicit evidence: InfField[A]): CFastDBM[M, SM, Closed, A] =
        dbm match {
          case BottomFast(n) => BottomFast(n)
          case CFast(m: FastDBM[M, SM, A]) => CFast(m)
          case NCFast(m: FastDBM[M, SM, A]) => m.incrementalClosure(v)
        }

      def strongClosure[S <: DBMState, A](dbm: CFastDBM[M, SM, S, A])
                          (implicit evidence: InfField[A]): CFastDBM[M, SM, Closed, A] =
        dbm match {
          case BottomFast(n) => BottomFast(n)
          case CFast(m: FastDBM[M, SM, A]) => CFast(m)
          case NCFast(m: FastDBM[M, SM, A]) => m.strongClosure
        }

      def forget[A](vi: VarIndex)(m: CFastDBM[M, SM, Closed, A])
        (implicit f: InfField[A]): CFastDBM[M, SM, Closed, A] =
        // TODO: this is eccessively approximate for decomposed matrices.
        // Forget operator for those should remove the forgot variable from the
        // indep. component including it (if any).
        mapFastDBM[M, SM, Closed, A](mapInnerMatrix(mev.ds.forget(vi)))(m)

      def nOfVars[S <: DBMState, A](m: CFastDBM[M, SM, S, A]): VarCount =
        FastDBMUtils.nOfVars(m)

      def get[S <: DBMState, A](i: Int, j: Int)(m: CFastDBM[M, SM, S, A])
                               (implicit ifield: InfField[A]): Option[A] =
        cfastInnerMatrix(m).map((inner) => mev.ds.get(i, j)(inner))

      def dbmIntersection[A, R <: DBMState, S <: DBMState]
        (m1: CFastDBM[M, SM, R, A], m2: CFastDBM[M, SM, S, A])
        (implicit ifield: InfField[A]): ExistsM[A] = {

        (m1, m2) match {
          case (BottomFast(nOfVars), _) =>
            packEx(BottomFast(nOfVars))
          case (_, BottomFast(nOfVars)) =>
            packEx(BottomFast(nOfVars))
          case (CFast(dbm1), CFast(dbm2)) =>
            packEx(NCFast(dbm1.intersection(dbm2)))
          case (CFast(dbm1), NCFast(dbm2)) =>
            packEx(NCFast(dbm1.intersection(dbm2)))
          case (NCFast(dbm1), CFast(dbm2)) =>
            packEx(NCFast(dbm1.intersection(dbm2)))
          case (NCFast(dbm1), NCFast(dbm2)) =>
            packEx(NCFast(dbm1.intersection(dbm2)))
        }
      }

      def topDBM[A](nOfVars: VarCount)(implicit ifield: InfField[A]):
          CFastDBM[M, SM, Closed, A] = {
        val top = mev.dec.pure(nOfVars, ifield.infinity)
        val cTop = allIndices(varCountToDim(nOfVars))
          .foldLeft(top)((m, i) => mev.ds.update(i, i, ifield.zero)(m))
        CFast(FullDBM(cTop, mev))
      }

      def bottomDBM[A](nOfVars: VarCount)(implicit ifield: InfField[A]):
          CFastDBM[M, SM, Closed, A] = BottomFast(nOfVars)

      def fromFun[A](d: Dimension, f: ((Int, Int) => A))
        (implicit ifield: InfField[A]): CFastDBM[M, SM, NonClosed, A] =
        NCFast(FullDBM(mev.ds.update(f)(
          mev.dec.pure(dimToVarCount(d), ifield.infinity)), mev))

      def flipVar[S <: DBMState, A](vi: VarIndex)(m: CFastDBM[M, SM, S, A])
                                   (implicit ifield: InfField[A]): CFastDBM[M, SM, S, A] =
        mapFastDBM[M, SM, S, A](mapInnerMatrix(mev.ds.flipVar(vi)))(m)

      def dbmUnion[A]
        (m1: CFastDBM[M, SM, Closed, A], m2: CFastDBM[M, SM, Closed, A])
        (implicit ifield: InfField[A]): CFastDBM[M, SM, Closed, A] = {

        (m1, m2) match {
          case (BottomFast(nOfVars), _) => BottomFast(nOfVars)
          case (_, BottomFast(nOfVars)) => BottomFast(nOfVars)
          case (CFast(dbm1), CFast(dbm2)) => CFast(dbm1.union(dbm2))
        }
      }

      def addScalarOnVar[S <: DBMState, A]
        (vi: VarIndex, const: A)(m: CFastDBM[M, SM, S, A])
        (implicit ifield: InfField[A]): CFastDBM[M, SM, S, A] =
        mapFastDBM[M, SM, S, A](
          mapInnerMatrix(mev.ds.addScalarOnVar(vi, const)))(m)

      def isBottomDBM[A, S <: DBMState](m: CFastDBM[M, SM, S, A])
                                      (implicit ifield: InfField[A]): Boolean =
        m match { case BottomFast(_) => true ; case _ => false }

      def widening[A, R <: DBMState, S <: DBMState]
        (dbm1: CFastDBM[M, SM, R, A], dbm2: CFastDBM[M, SM, S, A])
        (implicit ifield: InfField[A]): ExistsM[A] = {

        (dbm1, dbm2) match {
          case (BottomFast(nOfVars), _) => packEx(BottomFast(nOfVars))
          case (_, BottomFast(nOfVars)) => packEx(BottomFast(nOfVars))
          case (CFast(m1), CFast(m2))   => packEx(NCFast(m1.widening(m2)))
          case (CFast(m1), NCFast(m2))  => packEx(NCFast(m1.widening(m2)))
          case (NCFast(m1), CFast(m2))  => packEx(NCFast(m1.widening(m2)))
          case (NCFast(m1), NCFast(m2)) => packEx(NCFast(m1.widening(m2)))
        }
      }

      def narrowing[A, R <: DBMState, S <: DBMState]
        (dbm1: CFastDBM[M, SM, R, A], dbm2: CFastDBM[M, SM, S, A])
        (implicit ifield: InfField[A]): ExistsM[A] = {

        val nOfVars = FastDBMUtils.nOfVars(dbm1)
        (cfastInnerMatrix(dbm1), cfastInnerMatrix(dbm2)) match {
          case (None, _) => packEx(BottomFast(nOfVars))
          case (_, None) => packEx(BottomFast(nOfVars))
          // TODO: why not component-wise operator for decomposed matrices?
          case (Some(m1), Some(m2)) =>
            val newMat = mev.ds.narrowing(m1, m2)
            packEx(NCFast(FullDBM(newMat, mev)))
        }
      }

      def isTopDBM[A, S <: DBMState](dbm: CFastDBM[M,SM, S,A])
        (implicit ifield: InfField[A]): Boolean =
        cfastInnerMatrix(dbm) match {
          case Some(m) =>
            grid(varCountToDim(mev.ds.nOfVars(m)))
              .forall({ case (i, j) => mev.ds.get(i, j)(m) == ifield.infinity })
          case None => false
        }

      def addVariable[S <: DBMState, A](dbm: CFastDBM[M,SM, S,A])
          (implicit ifield: InfField[A]): CFastDBM[M,SM, S,A] = dbm match {
            case BottomFast(n) => BottomFast(addOne(n))
            case _ =>
              mapFastDBM[M, SM, S, A](
                mapInnerMatrix[M, SM, A](mev.dec.addVariable))(dbm)
          }

      def decideState[S <: DBMState, A](dbm: CFastDBM[M,SM, S,A]):
          DBMIxed[({ type T[W,B] = CFastDBM[M, SM, W, B]})#T, A] =  dbm match {
          case m @ BottomFast(_) =>
            CIxed[({ type T[W,B] = CFastDBM[M, SM, W, B]})#T, A](m)
          case m @ CFast(_) =>
            CIxed[({ type T[W,B] = CFastDBM[M, SM, W, B]})#T, A](m)
          case m @ NCFast(_) =>
            NCIxed[({ type T[W,B] = CFastDBM[M, SM, W, B]})#T, A](m)
        }

      def deleteVariable[S <: DBMState, A](v: VarIndex)(dbm: CFastDBM[M,SM,S,A])
        (implicit ifield: InfField[A]): CFastDBM[M,SM,S,A] = {
        def f(w: VarIndex): Option[VarIndex] =
          if (w == v) None else
            if (w < v) Some(w) else
              Some(VarIndex(w.i - 1))
        mapVariables(f)(dbm)
      }

      def mapVariables[S <: DBMState, A]
        (f: VarIndex => Option[VarIndex])(dbm: CFastDBM[M,SM,S,A])
        (implicit ifield: InfField[A]): CFastDBM[M,SM,S,A] = dbm match {
        case BottomFast(n) =>
          val newN = allVars(n).count(f(_).isDefined)
          BottomFast(VarCount(newN))
        case _ =>
          mapFastDBM[M, SM, S, A](fast =>
            mapInnerMatrix[M, SM, A](mev.dec.mapVariables(f))(fast.toFull))(dbm)
      }

      def compare[A](x: ExistsM[A], y: ExistsM[A])
                    (implicit evidence: InfField[A]): Option[Ordering] = {
        val dbm1: Option[M[A]] = cfastInnerMatrix(x.elem)
        val dbm2: Option[M[A]] = cfastInnerMatrix(y.elem)
        (dbm1, dbm2) match {
          case (None, None) => Some(EQ)
          case (Some(_), None) => Some(GT)
          case (None, Some(_)) => Some(LT)
          case (Some(m1), Some(m2)) => mev.ds.compare(m1, m2)
        }
      }
    }

  def packEx[M[_], SM[_], S <: DBMState, A](fastDBM: CFastDBM[M, SM, S, A]):
      ExistsDBM[({ type T[W] = CFastDBM[M, SM, W, A]})#T] =
    MkEx[S, ({ type T[S] = CFastDBM[M, SM, S, A]})#T](fastDBM)

  def cfastInnerMatrix[M[_], SM[_], S <: DBMState, A]
    (cfdbm: CFastDBM[M, SM, S, A])
    (implicit mev: MEvidence[M, SM], ifield: InfField[A]): Option[M[A]] =
    cfdbm match {
      case CFast(m: FastDBM[M, SM, A]) => Some(FastDBMUtils.fastInnerMatrix(m))
      case NCFast(m: FastDBM[M, SM, A]) => Some(FastDBMUtils.fastInnerMatrix(m))
      case BottomFast(_) => None
    }

  def mapFastDBM[M[_], SM[_], S <: DBMState, A](f: FastDBM[M, SM, A] => FastDBM[M, SM, A])
    (cfdbm: CFastDBM[M, SM, S, A])(implicit mev: MEvidence[M, SM], ifield: InfField[A]):
      CFastDBM[M, SM, S, A] =
    cfdbm match {
      case CFast(m: FastDBM[M, SM, A]) => CFast(f(m))
      case NCFast(m: FastDBM[M, SM, A]) => NCFast(f(m))
      case BottomFast(n) => BottomFast(n)
    }

  def mapInnerMatrix[M[_], SM[_], A](f: M[A] => M[A])(dbm: FastDBM[M, SM, A])
    (implicit mev: MEvidence[M, SM], ifield: InfField[A]): FastDBM[M, SM, A] =
    dbm match {
      case FullDBM(m, _) => FullDBM(f(m), mev)
      case DecomposedDBM(m, comps, _) => DecomposedDBM(f(m), comps, mev)
    }

}
