package it.unich.jandom.domains.numerical.octagon.fast

import breeze.numerics.pow
import it.unich.jandom.domains.numerical.octagon._
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

object CFDBMInstance {
  def instance[M[_]](implicit ds: DenseSparseDBM[M]) =
    new DifferenceBoundMatrix[({ type T[S, A] = CFastDBM[M, S, A] })#T] {

      type PosetConstraint[A] = InfField[A]

      def update[S <: DBMState, A](f: (Int, Int) => A)(m: CFastDBM[M, S, A])
                                  (implicit ifield: InfField[A]): ExistsM[A] =
        Utils.liftFromInner(ds.update(f))(m)

      def incrementalClosure[S <: DBMState, A](v: VarIndex)
                               (dbm: CFastDBM[M, S, A])
        (implicit evidence: InfField[A]): CFastDBM[M, Closed, A] =
        dbm match {
          case BottomFast(n) => BottomFast(n)
          case TopFast(n) => TopFast(n)
          case CFast(m) => CFast(m)
          case NCFast(m) => m.incrementalClosure(v)
        }

      def strongClosure[S <: DBMState, A](dbm: CFastDBM[M, S, A])
                          (implicit evidence: InfField[A]): CFastDBM[M, Closed, A] =
        dbm match {
          case BottomFast(n) => BottomFast(n)
          case TopFast(n) => TopFast(n)
          case CFast(m) => CFast(m)
          case NCFast(m) => m.strongClosure
        }

      def forget[S <: DBMState, A](vi: VarIndex)(m: CFastDBM[M, S, A])
                                  (implicit ifield: InfField[A]): CFastDBM[M, S, A] =
        Utils.mapFastDBM[M, S, A](fast =>
          Utils.mapInnerMatrix[M, A](inner =>
            ds.forget(vi)(inner)
          )(fast)
        )(m)

      def nOfVars[S <: DBMState, A](m: CFastDBM[M, S, A]): VarCount = Utils.nOfVars(m)

      def get[S <: DBMState, A](i: Int, j: Int)(m: CFastDBM[M, S, A])
                               (implicit ifield: InfField[A]): Option[A] =
        Utils.cfastInnerMatrix(m).flatMap((inner) => ds.get(i, j)(inner))

      def dbmIntersection[A, R <: DBMState, S <: DBMState]
        (m1: CFastDBM[M, R, A], m2: CFastDBM[M, S, A])
        (implicit ifield: InfField[A]): ExistsM[A] = {

        def aux(m1: FastDBM[M, A], m2: FastDBM[M, A]): FastDBM[M, A] =
          (m1, m2) match {
            case (FullDBM(dbm1, dsdbm), FullDBM(dbm2, _)) =>
              FullDBM(dsdbm.dbmUnion(dbm1, dbm2), dsdbm)
            case (DecomposedDBM(dbm1, comps1, dsdbm), DecomposedDBM(dbm2, comps2, _)) =>
              // compute sets of initialized variables, i.e., variables
              // that can be different from infinity
              val vars1 = comps1.foldRight(Seq[VarIndex]())(_ ++ _).toSet
              val vars2 = comps2.foldRight(Seq[VarIndex]())(_ ++ _).toSet

              // create a new component with the variables that can be
              // different from infinity in the intersection
              val vars = vars1 union vars2
              val component = vars.toSeq

              // create new submatrices with the same components
              val sub1 = dsdbm.extract(component)(dbm1)
              val sub2 = dsdbm.extract(component)(dbm2)

              val newMat = dsdbm.dbmIntersection(sub1, sub2)
              DecomposedDBM(newMat, Seq(component), dsdbm)
            case (dbm1 @ DecomposedDBM(_, _, _), dbm2 @ FullDBM(_, _)) =>
              aux(dbm1.toFull, dbm2)
            case (dbm1 @ FullDBM(_, _), dbm2 @ DecomposedDBM(_, _, _)) =>
              aux(dbm1, dbm2.toFull)
          }

        (m1, m2) match {
          case (BottomFast(nOfVars), _) =>
            Utils.packEx(BottomFast(nOfVars))
          case (_, BottomFast(nOfVars)) =>
            Utils.packEx(BottomFast(nOfVars))
          case (TopFast(_), other) =>
            Utils.packEx(other)
          case (other, TopFast(_)) =>
            Utils.packEx(other)
          case (CFast(dbm1), CFast(dbm2)) =>
            Utils.packEx(NCFast(aux(dbm1, dbm2)))
          case (CFast(dbm1), NCFast(dbm2)) =>
            Utils.packEx(NCFast(aux(dbm1, dbm2)))
          case (NCFast(dbm1), CFast(dbm2)) =>
            Utils.packEx(NCFast(aux(dbm1, dbm2)))
          case (NCFast(dbm1), NCFast(dbm2)) =>
            Utils.packEx(NCFast(aux(dbm1, dbm2)))
        }
      }

      def topDBM[A](nOfVars: VarCount)(implicit ifield: InfField[A]): CFastDBM[M, Closed, A] =
        TopFast(nOfVars)

      def bottomDBM[A](nOfVars: VarCount)(implicit ifield: InfField[A]): CFastDBM[M, Closed, A] =
        BottomFast(nOfVars)

      def fromFun[A](d: Dimension, f: ((Int, Int) => A))(implicit ifield: InfField[A]): CFastDBM[M, Closed, A] =
        CFast(FullDBM(ds.update(f)(ds.pure(halvedDimension(d), ifield.infinity)), ds))

      def flipVar[S <: DBMState, A](vi: VarIndex)(m: CFastDBM[M, S, A])
                                   (implicit ifield: InfField[A]): CFastDBM[M, S, A] =
        Utils.mapFastDBM[M, S, A](fast =>
          Utils.mapInnerMatrix[M, A](inner =>
            ds.flipVar(vi)(inner)
          )(fast)
        )(m)

      def dbmUnion[S <: DBMState, A](m1: CFastDBM[M, S, A], m2: CFastDBM[M, S, A])
                                    (implicit ifield: InfField[A]): CFastDBM[M, S, A] = {

        def aux(m1: FastDBM[M, A], m2: FastDBM[M, A]): FastDBM[M, A] =
          (m1, m2) match {
            case (FullDBM(dbm1, dsdbm), FullDBM(dbm2, _)) =>
              FullDBM(dsdbm.dbmUnion(dbm1, dbm2), dsdbm)
            case (DecomposedDBM(dbm1, comps1, dsdbm), DecomposedDBM(dbm2, comps2, _)) =>
              val vars1 = comps1.foldRight(Seq[VarIndex]())(_ ++ _).toSet
              val vars2 = comps2.foldRight(Seq[VarIndex]())(_ ++ _).toSet
              val vars = vars1 intersect vars2
              val newComps = comps1.map(_.filter(vars.contains(_)))
              val matrices = newComps.map(c => {
                  val sub1 = dsdbm.extract(c)(dbm1)
                  val sub2 = dsdbm.extract(c)(dbm2)
                  dsdbm.dbmUnion(sub1, sub2)
                })
              val newMat = matrices.foldLeft(dbm1)((mat, subMat) =>
                  dsdbm.pour(subMat)(mat)
                )
              DecomposedDBM(newMat, newComps, dsdbm)
            case (dbm1 @ DecomposedDBM(_, _, _), dbm2 @ FullDBM(_, _)) =>
              aux(dbm1.toFull, dbm2)
            case (dbm1 @ FullDBM(_, _), dbm2 @ DecomposedDBM(_, _, _)) =>
              aux(dbm1, dbm2.toFull)
          }

        (m1, m2) match {
          case (BottomFast(nOfVars), _) =>
            BottomFast(nOfVars)
          case (_, BottomFast(nOfVars)) =>
            BottomFast(nOfVars)
          case (TopFast(nOfVars), _) =>
            TopFast(nOfVars)
          case (_, TopFast(nOfVars)) =>
            TopFast(nOfVars)
          case (CFast(dbm1), CFast(dbm2)) =>
            CFast(aux(dbm1, dbm2))
          case (CFast(dbm1), NCFast(dbm2)) =>
            // WTF doesn't work TODO
            // NCFast(aux(dbm1, dbm2))
            dbmUnion(m2, m1)
          case (NCFast(dbm1), CFast(dbm2)) =>
            NCFast(aux(dbm1, dbm2))
          case (NCFast(dbm1), NCFast(dbm2)) =>
            NCFast(aux(dbm1, dbm2))
        }
      }

      def addScalarOnVar[S <: DBMState, A](vi: VarIndex, const: A)
                                          (m: CFastDBM[M, S, A])
                                          (implicit ifield: InfField[A]): CFastDBM[M, S, A] =
        Utils.mapFastDBM[M, S, A](fast =>
          Utils.mapInnerMatrix[M, A](inner =>
            ds.addScalarOnVar(vi, const)(inner)
          )(fast)
        )(m)

      def isBottomDBM[A, S <: DBMState](m: CFastDBM[M, S, A])
                                      (implicit ifield: InfField[A]): Boolean =
        Utils.cfastInnerMatrix(m).isEmpty

      def widening[A, R <: DBMState, S <: DBMState]
        (dbm1: CFastDBM[M, R, A], dbm2: CFastDBM[M, S, A])
        (implicit ifield: InfField[A]): ExistsM[A] = {

        def aux(m1: FastDBM[M, A], m2: FastDBM[M, A]): FastDBM[M, A] =
          (m1, m2) match {
            case (FullDBM(dbm1, dsdbm), FullDBM(dbm2, _)) =>
              FullDBM(dsdbm.widening(dbm1, dbm2), dsdbm)
            case (DecomposedDBM(dbm1, comps1, dsdbm), DecomposedDBM(dbm2, comps2, _)) =>
              val vars1 = comps1.foldRight(Seq[VarIndex]())(_ ++ _).toSet
              val vars2 = comps2.foldRight(Seq[VarIndex]())(_ ++ _).toSet
              val vars = vars1 intersect vars2
              val newComps = comps1.map(_.filter(vars.contains(_)))
              val matrices = newComps.map(c => {
                  val sub1 = dsdbm.extract(c)(dbm1)
                  val sub2 = dsdbm.extract(c)(dbm2)
                  dsdbm.widening(sub1, sub2)
                })
              val newMat = matrices.foldLeft(dbm1)((mat, subMat) =>
                  dsdbm.pour(subMat)(mat)
                )
              DecomposedDBM(newMat, newComps, dsdbm)
            case (dbm1 @ DecomposedDBM(_, _, _), dbm2 @ FullDBM(_, _)) =>
              aux(dbm1.toFull, dbm2)
            case (dbm1 @ FullDBM(_, _), dbm2 @ DecomposedDBM(_, _, _)) =>
              aux(dbm1, dbm2.toFull)
          }

        (dbm1, dbm2) match {
          case (BottomFast(nOfVars), _) =>
            Utils.packEx(BottomFast(nOfVars))
          case (_, BottomFast(nOfVars)) =>
            Utils.packEx(BottomFast(nOfVars))
          case (TopFast(nOfVars), _) =>
            Utils.packEx(TopFast(nOfVars))
          case (_, TopFast(nOfVars)) =>
            Utils.packEx(TopFast(nOfVars))
          case (CFast(m1), CFast(m2)) =>
            Utils.packEx(NCFast(aux(m1, m2)))
          case (CFast(m1), NCFast(m2)) =>
            Utils.packEx(NCFast(aux(m1, m2)))
          case (NCFast(m1), CFast(m2)) =>
            Utils.packEx(NCFast(aux(m1, m2)))
          case (NCFast(m1), NCFast(m2)) =>
            Utils.packEx(NCFast(aux(m1, m2)))
        }
      }

      def narrowing[A, R <: DBMState, S <: DBMState]
        (dbm1: CFastDBM[M, R, A], dbm2: CFastDBM[M, S, A])
        (implicit ifield: InfField[A]): ExistsM[A] = {
        val nOfVars = Utils.nOfVars(dbm1)
        (Utils.cfastInnerMatrix(dbm1), Utils.cfastInnerMatrix(dbm2)) match {
          case (None, _) => Utils.packEx(BottomFast(nOfVars))
          case (_, None) => Utils.packEx(BottomFast(nOfVars))
          case (Some(m1), Some(m2)) =>
            val newMat = ds.narrowing(m1, m2)
            Utils.packEx(NCFast(FullDBM(newMat, ds)))
        }
      }

      def isTopDBM[A, S <: DBMState](dbm: CFastDBM[M,S,A])(implicit ifield: InfField[A]): Boolean =
        dbm match {
          case TopFast(_) => true
          case _ => false
        }

      def addVariable[S <: DBMState, A](dbm: CFastDBM[M,S,A])
          (implicit ifield: InfField[A]): CFastDBM[M,S,A] = dbm match {
            case BottomFast(n) => BottomFast(addOne(n))
            case TopFast(n) => TopFast(addOne(n))
            case m =>
              Utils.mapFastDBM[M, S, A](fast =>
                Utils.mapInnerMatrix[M, A](inner =>
                  ds.addVariable(inner)
                )(fast)
              )(dbm)
          }

      def decideState[S <: DBMState, A](dbm: CFastDBM[M,S,A]):
          DBMIxed[({ type T[W,B] = CFastDBM[M, W, B]})#T, A] =  dbm match {
          case m @ BottomFast(_) =>
            CIxed[({ type T[W,B] = CFastDBM[M, W, B]})#T, A](m)
          case m @ TopFast(_) =>
            CIxed[({ type T[W,B] = CFastDBM[M, W, B]})#T, A](m)
          case m @ CFast(_) =>
            CIxed[({ type T[W,B] = CFastDBM[M, W, B]})#T, A](m)
          case m @ NCFast(_) =>
            NCIxed[({ type T[W,B] = CFastDBM[M, W, B]})#T, A](m)
        }

      def deleteVariable[S <: DBMState, A](v: VarIndex)(dbm: CFastDBM[M,S,A])
          (implicit ifield: InfField[A]): CFastDBM[M,S,A] = dbm match {
            case BottomFast(n) => BottomFast(subOne(n))
            case TopFast(n) => TopFast(subOne(n))
            case m =>
              Utils.mapFastDBM[M, S, A](fast =>
                Utils.mapInnerMatrix[M, A](inner =>
                  ds.deleteVariable(inner)
                )(fast)
              )(dbm)
          }

      def mapVariables[S <: DBMState, A](f: VarIndex => Option[VarIndex])
          (dbm: CFastDBM[M,S,A])(implicit ifield: InfField[A]): CFastDBM[M,S,A] = dbm match {
        case BottomFast(n) =>
          val newN = allVars(n).count(f(_).isDefined)
          BottomFast(VarCount(newN))
        case TopFast(n) =>
          val newN = allVars(n).count(f(_).isDefined)
          BottomFast(VarCount(newN))
        case m =>
          Utils.mapFastDBM[M, S, A](fast =>
            Utils.mapInnerMatrix[M, A](inner =>
              ds.mapVariables(f)(inner)
            )(fast)
          )(dbm)
      }

      def compare[A](x: ExistsM[A], y: ExistsM[A])
                    (implicit evidence: InfField[A]): Option[Ordering] = {
        val dbm1: Option[M[A]] = Utils.cfastInnerMatrix(x.elem)
        val dbm2: Option[M[A]] = Utils.cfastInnerMatrix(y.elem)
        (dbm1, dbm2) match {
          case (None, None) => Some(EQ)
          case (Some(_), None) => Some(GT)
          case (None, Some(_)) => Some(LT)
          case (Some(m1), Some(m2)) => ds.compare(m1, m2)
        }
      }
    }
}

// ADT of "closable" DBMs in their fast implementation from Vechev et al.
// They are "closable" in the sense that they augment the ADT of fast DBMs with
// the type-level capability of being indexed by their strong closure state.
sealed trait CFastDBM[M[_], _, A]
// Constructor of *closed* fast DBMs.
case class CFast[M[_], A](m: FastDBM[M, A]) extends CFastDBM[M, Closed, A]
// Constructor of *non-closed* fast DBMs.
case class NCFast[M[_], A](m: FastDBM[M, A]) extends CFastDBM[M, NonClosed, A]
case class TopFast[M[_], A](nOfVars: VarCount) extends CFastDBM[M, Closed, A]
case class BottomFast[M[_], A](nOfVars: VarCount) extends CFastDBM[M, Closed, A]

object Utils {

  def packEx[M[_], S <: DBMState, A](fastDBM: CFastDBM[M, S, A])
  : ExistsDBM[({ type T[W] = CFastDBM[M, W, A]})#T] =
    MkEx[S, ({ type T[S] = CFastDBM[M, S, A]})#T](fastDBM)

  def nOfVars[M[_], S, A](dbm: CFastDBM[M, S, A])
                         (implicit ds: DenseSparseDBM[M]): VarCount =
    dbm match {
      case CFast(m: FastDBM[M, A]) => ds.nOfVars(fastInnerMatrix(m))
      case NCFast(m: FastDBM[M, A]) => ds.nOfVars(fastInnerMatrix(m))
      case TopFast(n) => n
      case BottomFast(n) => n
    }

  def fastInnerMatrix[M[_], S <: DBMState, A](fdbm: FastDBM[M, A]): M[A] =
    fdbm match {
      case FullDBM(m: M[A], _) => m
      case DecomposedDBM(m: M[A], _, _) => m
    }

  def cfastInnerMatrix[M[_], S <: DBMState, A]
    (cfdbm: CFastDBM[M, S, A])(implicit ds: DenseSparseDBM[M], ifield: InfField[A]): Option[M[A]] = {

    cfdbm match {
      case CFast(m: FastDBM[M, A]) => Some(fastInnerMatrix(m))
      case NCFast(m: FastDBM[M, A]) => Some(fastInnerMatrix(m))
      case TopFast(n) => Some(ds.pure(n, ifield.infinity))
      case BottomFast(_) => None
    }
  }

  // Applying an arbitrary map destroys every knowledge about relations between
  // variables, so we just return DBM with unknown closure state (because it
  // could be either non-closed or bottom, which is closed).
  // WARNING: input and output matrices of f as assumed to be of the same
  // dimension!
  def liftFromInner[M[_], S <: DBMState, A](f : M[A] => M[A])
                                           (dbm: CFastDBM[M, S, A])
                                           (implicit ds: DenseSparseDBM[M],
                                            ifield: InfField[A])
  : ExistsDBM[({ type T[W] = CFastDBM[M, W, A]})#T] =
    Utils.cfastInnerMatrix(dbm).map((inner) =>
      NCFast(FullDBM(f(inner), ds))) match {
      case Some(cm) => packEx(cm)
      case None => packEx(BottomFast(nOfVars(dbm)))
    }

  def mapFastDBM[M[_], S <: DBMState, A](f: FastDBM[M, A] => FastDBM[M, A])
    (cfdbm: CFastDBM[M, S, A])(implicit ds: DenseSparseDBM[M], ifield: InfField[A]): CFastDBM[M, S, A] = {

    cfdbm match {
      case CFast(m: FastDBM[M, A]) => CFast(f(m))
      case NCFast(m: FastDBM[M, A]) => NCFast(f(m))
      case TopFast(n) => CFast(f(FullDBM(ds.pure(n, ifield.infinity), ds)))
      case BottomFast(n) => BottomFast(n)
    }
  }

  def mapInnerMatrix[M[_], A](f: M[A] => M[A])(dbm: FastDBM[M, A])
                             (implicit ds: DenseSparseDBM[M], ifield: InfField[A]): FastDBM[M, A] = {
    dbm match {
      case FullDBM(m, _) => FullDBM(f(m), ds)
      case DecomposedDBM(m, comps, _) => DecomposedDBM(f(m), comps, ds)
    }
  }

}

object FastDbmUtils {
  val sparseThreshold = 0.5

  def nuffDecomposed(is: List[List[VarIndex]]): Boolean = is.size > 1
  def nuffSparse(d: VarCount, is: NNI): Boolean =
    1.0 - (is.nni / (2 * d.count * d.count + 2 * d.count)) >= sparseThreshold

  def calculateComponents[M[_], A](dbm: FastDBM[M, A])
                                   (implicit e: DenseSparseDBM[M],
                                    ifield: InfField[A]): List[List[VarIndex]] = {
    val innerMatrix = Utils.fastInnerMatrix(dbm)
    def related(vi: VarIndex, vj: VarIndex): Boolean = {
      import VarIndexOps._
      Set(
        (varPlus(vi), varPlus(vj)),
        (varPlus(vi), varMinus(vj)),
        (varMinus(vi), varPlus(vj)),
        (varMinus(vi), varMinus(vj))
      ).filter({ case (i, j) =>
        i != j
      }).exists({ case (i, j) =>
        e.get(i, j)(innerMatrix) match {
          case Some(v) => ifield.!=(v, ifield.infinity)
          case None    => false
        }
      })
    }

    import CountOps._

    val nOfVars = e.nOfVars(innerMatrix)
    val rels = for (vi <- allVars(nOfVars);
                    vj <- allVars(nOfVars);
                    if related(vi, vj))
                  yield (vi, vj)

    val indComps =
      rels.foldLeft(Set[Set[VarIndex]]())({ case (comps, (vi, vj)) =>
        val compI = comps.find(_.contains(vi)).getOrElse(Set())
        val compJ = comps.find(_.contains(vj)).getOrElse(Set())
        val newComp = Set(vi, vj) ++ compI ++ compJ
        comps.filter(c => !c.contains(vi) && !c.contains(vj)) + newComp
      })

    indComps.map(_.toList).toList
  }

}

sealed trait FastDBM[M[_], A] {
  // Skeleton of strong closure for fast matrices.
  def strongClosure(implicit rdbm: DenseSparseDBM[M], ifield: InfField[A])
  : CFastDBM[M, Closed, A] = {

    val dbm = Utils.fastInnerMatrix(this)
    val indepComponents: List[List[VarIndex]] =
          FastDbmUtils.calculateComponents(this)

    val submatrices = indepComponents.map(seq => rdbm.extract(seq)(dbm))
    Applicative[Option].sequence(
      submatrices.map(m => rdbm.strongClosure(m)).toList) match {

      case Some(closedSubs) => {
        val newMatrix = closedSubs.foldRight(dbm)({
            case (sub, full) => rdbm.pour(sub)(full)
          })

        if (FastDbmUtils.nuffDecomposed(indepComponents))
          CFast(DecomposedDBM(newMatrix, indepComponents, rdbm))
        else
          // might be that we have to compute the index of non-infinite terms
          // during the strong closure above, and pass them to the dbm
          // constructor.
          CFast(FullDBM(newMatrix, rdbm))
      }
      case None => BottomFast(VarCount(rdbm.varIndices(dbm).length))
    }
  }

    def incrementalClosure(v: VarIndex)(implicit ifield: InfField[A])
    : CFastDBM[M, Closed, A]
}

// Full DBMs are fast DBMs that are not decomposed, i.e., they can be either
// dense or sparse.
// Sparsity details, including when to switch between dense and sparse
// representation, is supposed to be handled by the specific implementation of
// the the DenseSparse trait/typeclass.
// An even better thing to do (time permitting) could be to use a suitable
// abstract trait of DBMs that does not talk about sparsity at all (who cares
// if the implementations use a dense/sparse representation anyway, as long as
// they provide the needed DBM-like operations?)
case class FullDBM[M[_], A](dbm: M[A], dsdbm: DenseSparseDBM[M]) extends FastDBM[M, A] {
  def incrementalClosure(v: VarIndex)(implicit ifield: InfField[A])
    : CFastDBM[M, Closed, A] = dsdbm.incrementalClosure(v)(dbm) match {
      case Some(m) => CFast(FullDBM(m, dsdbm))
      case None    => BottomFast(dsdbm.nOfVars(dbm))
    }
}

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
                                  indepComponents: Seq[Seq[VarIndex]],
                                  rdbm: DenseSparseDBM[M]) extends FastDBM[M, A] {

  def toFull: FullDBM[M, A] = FullDBM(completeDBM, rdbm)

  def incrementalClosure(v: VarIndex)(implicit ifield: InfField[A])
    : CFastDBM[M, Closed, A] = {
      val maybeComp = indepComponents.find(_.contains(v))
      maybeComp match {
        case Some(comp) =>
          val subMat = rdbm.extract(comp)(completeDBM)
          rdbm.incrementalClosure(v)(subMat) match {
            case Some(closed) =>
              val newMat = rdbm.pour(closed)(completeDBM)
              CFast(DecomposedDBM(newMat, indepComponents, rdbm))
            case None         =>
              BottomFast(rdbm.nOfVars(completeDBM))
          }
        case None       =>
          CFast(this)
      }
    }

}
