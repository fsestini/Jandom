package it.unich.jandom.domains.numerical.octagon
import it.unich.jandom.domains.numerical.octagon.variables.VarIndex
import it.unich.jandom.domains.numerical.octagon.variables.Dimension
import it.unich.jandom.domains.numerical.octagon.variables.VarCount
import it.unich.jandom.domains.numerical.octagon.variables.CountOps
import it.unich.jandom.domains.numerical.octagon.variables.VarIndexOps
import scala.language.higherKinds

// M-based raw DBM implementation
sealed trait DBM[M[A], S, A] {
  def liftFromInner(f: M[A] => M[A])(implicit ifield: InfField[A], me: Matrix[M]): DBM[M, S, A]
  def union(other: DBM[M, S, A])(implicit infField: InfField[A], me: Matrix[M]): DBM[M, S, A]
  def decideState: DBMIxed[({ type T[B,C] = ( DBM[M, B, C])})#T, A]
  val innerMatrix: Option[M[A]]
  def noOfVariables: VarCount
}

case class ClosedDBM[M[_], A](m: M[A])(implicit me: Matrix[M]) extends DBM[M, Closed, A] {
  type ThisDBM[B,C] = (DBM[M, B, C])
  def noOfVariables : VarCount = CountOps.dimToVarCount(me.dimension(m))
  override def liftFromInner(f: (M[A]) => M[A])
                            (implicit ifield: InfField[A], me: Matrix[M]): DBM[M, Closed, A] =
    ClosedDBM[M, A](f(m))

  override def union(other: DBM[M, Closed, A])
                    (implicit infField: InfField[A], me: Matrix[M]): DBM[M, Closed, A] = other match {
    case ClosedDBM(m2) => {
      require(noOfVariables == other.noOfVariables)
      ClosedDBM[M, A](me.combine(infField.max)(m, m2))
    }
    case BottomDBM(n2) => { require(noOfVariables == other.noOfVariables) ; this }
  }

  val innerMatrix: Option[M[A]] = Some(m)
  def decideState: DBMIxed[ThisDBM, A] = CIxed[ThisDBM,A](this)

  override def toString = "ClosedDBM[M]("+me.dimension(m)+") with M:\n" + m.toString
}

case class NonClosedDBM[M[_], A](m: M[A])(implicit me: Matrix[M]) extends DBM[M, NonClosed, A] {
  type ThisDBM[B,C] = (DBM[M, B, C])
  def noOfVariables : VarCount = CountOps.dimToVarCount(me.dimension(m))
  override def liftFromInner(f: (M[A]) => M[A])
                            (implicit ifield: InfField[A], me: Matrix[M]): DBM[M, NonClosed, A] =
    NonClosedDBM(f(m))
  def union(other: DBM[M, NonClosed, A])
           (implicit infField: InfField[A], me: Matrix[M]): DBM[M, NonClosed, A] = other match {
    case NonClosedDBM(m2) => {
      require(noOfVariables == other.noOfVariables)
      NonClosedDBM(me.combine(infField.max)(m, m2))
    }
  }

  val innerMatrix: Option[M[A]] = Some(m)

  def decideState: DBMIxed[ThisDBM, A] = NCIxed[ThisDBM, A](this)

  override def toString = "NonClosedDBM[M]("+me.dimension(m)+") with M:\n" + m.toString
}

case class BottomDBM[M[_], A](noOfVariables: VarCount) extends DBM[M, Closed, A] {
  type ThisDBM[B,C] = (DBM[M, B, C])
  override def liftFromInner(f: (M[A]) => M[A])
                            (implicit ifield: InfField[A], me: Matrix[M]): DBM[M, Closed, A] =
    BottomDBM[M, A](noOfVariables)

  def union(other: DBM[M, Closed, A])(implicit infField: InfField[A], me: Matrix[M]): DBM[M, Closed, A] =
    { require(noOfVariables == other.noOfVariables) ; other }

  val innerMatrix: Option[M[A]] = None

  def decideState: DBMIxed[ThisDBM, A] = {
    val dbm: DBM[M, Closed, A] = BottomDBM[M, A](noOfVariables)
    CIxed[ThisDBM, A](dbm)
  }

  override def toString = "BottomDBM[M]("+noOfVariables+")"
}

class DBMInstance[M[_]](implicit me: Matrix[M]) {
  type ThisDBM[B,C] = (DBM[M, B, C])

  implicit val funDBM: DifferenceBoundMatrix[ThisDBM] { type PosetConstraint[A] = InfField[A] } =
    new DifferenceBoundMatrix[ThisDBM] {
    def update[S <: DBMState, A](f: (Int, Int) => A)
                                (dbm: DBM[M, S, A])
                                (implicit ifield: InfField[A]): ExistsM[A] =
      dbm.innerMatrix match {
        case Some(inner) => mkExFun(NonClosedDBM[M,A](me.update(f)(inner)))
        case None => mkExFun(BottomDBM[M,A](dbm.noOfVariables))
      }

    def incrementalClosure[S <: DBMState, A](v: VarIndex)
                             (dbm: DBM[M, S, A])
      (implicit evidence: InfField[A]): DBM[M, Closed, A] =
      dbm.innerMatrix match {
        case Some(m) =>
          (new BagnaraStrongClosure[M,A]()).incrementalClosure(v)(m) match {
            case Some(closed: M[A]) => ClosedDBM[M, A](closed)
            case None => BottomDBM[M, A](dbm.noOfVariables)
          }
        case None => BottomDBM[M, A](dbm.noOfVariables)
      }

    def strongClosure[S <: DBMState, A](dbm: DBM[M, S, A])
                        (implicit evidence: InfField[A]): DBM[M, Closed, A] =
      dbm.innerMatrix match {
        case Some(m) =>
          (new BagnaraStrongClosure[M,A]()).strongClosure(m) match {
            case Some(closed: M[A]) => ClosedDBM[M,A](closed)
            case None => BottomDBM[M,A](dbm.noOfVariables)
          }
        case None => BottomDBM[M,A](dbm.noOfVariables)
      }

    def forget[S <: DBMState, A](vi: VarIndex)(m: DBM[M, S, A])
                                (implicit ifield: InfField[A]): DBM[M, S, A] = {
      m.liftFromInner((inner) => {
        val f: (Int, Int) => A = (i, j) => {
          if (i != VarIndexOps.varPlus(vi) && i != VarIndexOps.varMinus(vi) && j != VarIndexOps.varPlus(vi) && j != VarIndexOps.varMinus(vi))
            me.get(i, j)(inner)
          else if (i == j && (i == VarIndexOps.varPlus(vi) || i == VarIndexOps.varMinus(vi)))
            ifield.zero
          else ifield.infinity
        }
        me.update(f)(inner)
      })
    }

    def nOfVars[S <: DBMState, A](m: DBM[M, S, A]): VarCount = m.noOfVariables

    def get[S <: DBMState, A](i: Int, j: Int)(m: DBM[M, S, A])
                             (implicit ifield: InfField[A]): Option[A] =
      m.innerMatrix.map((mmm) => me.get(i, j)(mmm))

    def mkExFun[S <: DBMState, A](funDBM: DBM[M, S, A]): ExistsM[A] =
      MkEx[S, ({ type T[S] = DBM[M, S, A]})#T](funDBM)

    def dbmIntersection[A, S <: DBMState, T <: DBMState]
      (m1: DBM[M, S, A], m2: DBM[M, T, A])
      (implicit ifield: InfField[A]): ExistsM[A] = {
      require(m1.noOfVariables == m2.noOfVariables)
      val o = for {
        mm1 <- m1.innerMatrix
        mm2 <- m2.innerMatrix
      } yield NonClosedDBM[M,A](me.combine(ifield.min)(mm1, mm2))
      o match {
        case Some(matrix) => mkExFun(matrix)
        case None => mkExFun(BottomDBM[M,A](m1.noOfVariables))
      }
    }

    def topDBM[A](nOfVars: VarCount)(implicit ifield: InfField[A]): DBM[M, Closed, A] =
      ClosedDBM[M,A](me.make((i, j) => ifield.infinity, CountOps.varCountToDim(nOfVars)))

    def bottomDBM[A](nOfVars: VarCount)(implicit ifield: InfField[A]): DBM[M, Closed, A] = BottomDBM[M,A](nOfVars)
    def fromFun[A](d: Dimension, f: ((Int, Int) => A))(implicit ifield: InfField[A]): DBM[M, Closed, A] =
      strongClosure(NonClosedDBM[M,A](me.make(f, d)))
    def flipVar[S <: DBMState, A](vi: VarIndex)(dbm: DBM[M, S, A])
                                 (implicit ifield: InfField[A]): DBM[M, S, A] = {
      dbm.liftFromInner((inner) => {
        val f: (Int, Int) => A = (i, j) => {
          if (i == VarIndexOps.varPlus(vi) || i == VarIndexOps.varMinus(vi)) {
            if (j == VarIndexOps.varPlus(vi) || j == VarIndexOps.varMinus(vi))
              me.get(VarIndexOps.signed(i), VarIndexOps.signed(j))(inner)
            else
              me.get(VarIndexOps.signed(i), j)(inner)
          } else {
            if (j == VarIndexOps.varPlus(vi) || j == VarIndexOps.varMinus(vi))
              me.get(i, VarIndexOps.signed(j))(inner)
            else
              me.get(i, j)(inner)
          }
        }
        me.update(f)(inner)
      })
    }

    def dbmUnion[S <: DBMState, A](m1: DBM[M, S, A], m2: DBM[M, S, A])
                                  (implicit ifield: InfField[A]): DBM[M, S, A] = m1.union(m2)

    def addScalarOnVar[S <: DBMState, A](vi: VarIndex, const: A)
                                        (fundbm: DBM[M, S, A])
                                        (implicit ifield: InfField[A]): DBM[M, S, A] = {
      fundbm.liftFromInner((dbm) => {
        val f: (Int, Int) => A = (i, j) => {
          val g1 = (i == VarIndexOps.varPlus(vi) && j != VarIndexOps.varPlus(vi) && j != VarIndexOps.varMinus(vi)) ||
                   (j == VarIndexOps.varMinus(vi) && i != VarIndexOps.varPlus(vi) && i != VarIndexOps.varMinus(vi))
          val g2 = (i != VarIndexOps.varPlus(vi) && i != VarIndexOps.varMinus(vi) && j == VarIndexOps.varPlus(vi)) ||
                   (j != VarIndexOps.varPlus(vi) && j != VarIndexOps.varMinus(vi) && i == VarIndexOps.varMinus(vi))
          val g3 = i == VarIndexOps.varPlus(vi) && j == VarIndexOps.varMinus(vi)
          val g4 = i == VarIndexOps.varMinus(vi) && j == VarIndexOps.varPlus(vi)
          if (g1) ifield.-(me.get(i, j)(dbm), const) else
          if (g2) ifield.+(me.get(i, j)(dbm), const) else
          if (g3) ifield.-(me.get(i, j)(dbm), ifield.double(const)) else
          if (g4) ifield.+(me.get(i, j)(dbm), ifield.double(const)) else
            me.get(i, j)(dbm)
        }
        me.update(f)(dbm)
      })
    }

    def isBottomDBM[A, S <: DBMState](m: DBM[M, S, A])
                                     (implicit ifield: InfField[A]): Boolean =
      m match {
        case BottomDBM(_) => true
        case _ => false
      }

    def widening[A, S <: DBMState, T <: DBMState]
      (dbm1: DBM[M, S, A], dbm2: DBM[M, T, A])
      (implicit ifield: InfField[A]): ExistsM[A] = {

      require(dbm1.noOfVariables == dbm2.noOfVariables)


      (dbm1.innerMatrix, dbm2.innerMatrix) match {
        case (None, None)         => mkExFun(BottomDBM[M,A](dbm1.noOfVariables))
        case (Some(_), None)      => mkExFun(dbm1)
        case (None, Some(_))      => mkExFun(dbm2)
        case (Some(m1), Some(m2)) => mkExFun(NonClosedDBM[M,A](
                                          me.combine((mij: A, nij: A) =>
                                            ifield.compare(mij, nij) match {
                                              case GT => mij
                                              case _ => ifield.infinity
                                            }
                                          )(m1, m2)))
      }
    }

    def narrowing[A, S <: DBMState, T <: DBMState]
      (dbm1: DBM[M, S, A], dbm2: DBM[M, T, A])
      (implicit ifield: InfField[A]): ExistsM[A] = {

      require(dbm1.noOfVariables == dbm2.noOfVariables)

      val m = for {
        m1 <- dbm1.innerMatrix
        m2 <- dbm2.innerMatrix
      } yield me.combine((mij: A, nij: A) => if (mij == ifield.infinity) nij else mij)(m1, m2)
      m match {
        case Some(matrix) => mkExFun(NonClosedDBM[M,A](matrix))
        case None => mkExFun(BottomDBM[M,A](dbm1.noOfVariables))
      }
    }

    def decideState[S <: DBMState, A](dbm: DBM[M, S, A]): DBMIxed[ThisDBM, A] =
      dbm.decideState

    def isTopDBM[A, S <: DBMState](dbm: DBM[M, S, A])
                                  (implicit ifield: InfField[A]): Boolean =
      dbm.innerMatrix match {
        case Some(m) =>
          me.toList(m)
            .foldRight(true)((el, b) => el == ifield.infinity && b)
        case None => false
      }

    def addVariable[S <: DBMState, A](dbm: DBM[M, S, A])
                                     (implicit ifield: InfField[A]): DBM[M, S, A] =
      dbm.liftFromInner((m) => {
        me.make((i, j) => {
          if (CountOps.inDimension(i, j, CountOps.varCountToDim(dbm.noOfVariables)))
            me.get(i, j)(m)
          else ifield.infinity
        }, CountOps.varCountToDim(CountOps.addOne(dbm.noOfVariables)))
      })

    // Proved with pen and paper that shuffling variables around preserves
    // closure state.
    def mapVariables[S <: DBMState, A](f: (VarIndex) => Option[VarIndex])
                                      (dbm: DBM[M, S, A])
                                       (implicit ifield: InfField[A]): DBM[M, S, A] =
      dbm.liftFromInner((new VarMapping[M]).mapVariablesAux(f, dbm.noOfVariables))

    def deleteVariable[S <: DBMState, A](v: VarIndex)(dbm: DBM[M, S, A])
                                        (implicit ifield: InfField[A]): DBM[M, S, A] =
      mapVariables((x : VarIndex) =>
        if (x.i < v.i) Some(x) else
        if (x.i == v.i) None else Some(VarIndex(x.i - 1)))(dbm)

    type PosetConstraint[A] = InfField[A]

    def compare[A](x: ExistsDBM[({ type Q[S] = DBM[M, S, A]})#Q],
                   y: ExistsDBM[({ type Q[S] = DBM[M, S, A]})#Q])
                  (implicit evidence: InfField[A]): Option[Ordering] = {
      (x.elem.innerMatrix, y.elem.innerMatrix) match {
        case (None, None) => Some(EQ)
        case (Some(_), None) => Some(GT)
        case (None, Some(_)) => Some(LT)
        case (Some(m1), Some(m2)) => {
          val l: List[Ordering] = me.toList(me.combine(evidence.compare)(m1, m2))
          (l.forall(_ == EQ),
           l.forall(x => x == LT || x == EQ),
           l.forall(x => x == GT || x == EQ)) match {
            case (true, _, _) => Some(EQ)
            case (_, true, _) => Some(LT)
            case (_, _, true) => Some(GT)
            case _ => None
          }
        }
      }
    }

  }
}

// This is the simplified strong closure algorithm from
// Bagnara et al., Widening Operators for Weakly-Relational Numeric Abstractions.
//
// It is a classical Floyd-Warshall, followed by a single strong coherence step.
//
// for k = 0 to 2*n - 1
//   for i = 0 to 2*n - 1
//     for j = 0 to 2*n - 1
//       m(i,j) = min( m(i,j), m(i,k) + m(k,j) )
//
// for i = 0 to 2*n - 1
//   for j = 0 to 2*n - 1
//     m(i,j) = min( m(i,j), m(i, signed(i)) + m(signed(j), j) / 2)
class BagnaraStrongClosure[M[_], A](implicit val me: Matrix[M]) {

  private def signed(i: Int) = if (i % 2 == 0) i + 1 else i - 1

  def nullCheck[A](m: M[A])(implicit ifield: InfField[A]): Option[M[A]] = {
    val negative: Boolean = CountOps.allIndices(me.dimension(m)).exists((i) =>
      ifield.compare(me.get(i, i)(m), ifield.zero) == LT)
    if (negative) None else {
      val updater: (Int, Int) => A = (i, j) =>
        if (i == j) ifield.zero else me.get(i, j)(m)
      Some(me.update(updater)(m))
    }
  }

  def strengthen[A](dbm: M[A])(implicit ifield: InfField[A]): M[A] =
    CountOps.grid(me.dimension(dbm))
      .foldLeft(dbm)((x, pair) => pair match {
        case (i, j) => {
          val newVal =
            ifield.min(
              me.get(i, j)(x),
              ifield.half(ifield.+(me.get(i, signed(i))(x), me.get(signed(j), j)(x))))
          me.update(i, j, newVal)(x)
        }
      })

  def strongClosure[A](dbm: M[A])(implicit ifield: InfField[A]): Option[M[A]] = {
    val closed: M[A] = (for {
      k <- CountOps.allIndices(me.dimension(dbm))
      (i, j) <- CountOps.grid(me.dimension(dbm))
    } yield (k, i, j)).foldLeft(dbm)((x, triple) => triple match {
      case (k, i, j) => {
        val newVal =
          ifield.min(
            me.get(i,j)(x),
            ifield.+(me.get(i,k)(x), me.get(k,j)(x)))
        me.update(i, j, newVal)(x)
      }
    })

    nullCheck(strengthen(closed))
  }

  def incrementalClosure[A](vi: VarIndex)(dbm: M[A])
                           (implicit ifield: InfField[A]): Option[M[A]] = {
    val p: ((Int, Int, Int)) => Boolean = {
      case (k, i, j) =>
        k == vi.i || k == signed(vi.i) ||
          i == vi.i || i == signed(vi.i) ||
          j == vi.i || j == signed(vi.i)
    }

    val iclosed: M[A] =
      (for { k <- CountOps.allIndices(me.dimension(dbm)) ;
             (i, j) <- CountOps.grid(me.dimension(dbm)) } yield (k, i, j))
        .filter(p)
        .foldLeft(dbm)((x, triple) => triple match {
          case (k, i, j) => {
            val newVal =
              ifield.min(me.get(i,j)(x), ifield.+(me.get(i,k)(x), me.get(k,j)(x)))
            me.update(i, j, newVal)(x)
          }
        })

    nullCheck(strengthen(iclosed))
  }

}

class VarMapping[M[_]](implicit val me: Matrix[M]) {

  def varMapImageSize(f: VarIndex => Option[VarIndex], nOfVars: VarCount): VarCount =
    VarCount(
      CountOps.allVars(nOfVars).map((v) => f(v) match {
        case Some(_) => 1
        case None => 0
      }).sum)

  def mapVariablesAux[A](f: (VarIndex) => Option[VarIndex], nOfVars: VarCount)
                        (m: M[A]): M[A] = {

    def inverse(f: VarIndex => Option[VarIndex])(v: VarIndex): VarIndex = {
      val vars: Seq[VarIndex] = CountOps.allVars(nOfVars)
      val opts: Seq[Option[VarIndex]] = vars.map(x => f(x) match {
        case Some(vres) => if (v == vres) Some(x) else None
        case None => None
      })
      val opt: Option[VarIndex] = opts.flatten.headOption
      opt match {
        case Some(vres) => vres
        case None => throw new IllegalArgumentException("inverse out of bounds")
      }
    }

    def inverseIx(f: VarIndex => Option[VarIndex])(i: Int): Int =
      VarIndexOps.toIndexAndCoeff(i) match {
        case (v, VarIndexOps.Positive) => VarIndexOps.varPlus(inverse(f)(v))
        case (v, VarIndexOps.Negative) => VarIndexOps.varMinus(inverse(f)(v))
      }

    val newCount: VarCount = varMapImageSize(f, nOfVars)
    me.make((i, j) => {
      me.get(inverseIx(f)(i), inverseIx(f)(j))(m)
    }, CountOps.varCountToDim(newCount))
  }
}
