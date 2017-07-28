package it.unich.jandom.domains.numerical.octagon

import VarIndexOps._

// FunMatrix-based raw DBM implementation
sealed trait FunDBM[S, A] {
  def liftFromInner(f: FunMatrix[A] => FunMatrix[A])(implicit ifield: InfField[A]): FunDBM[S, A]
  def union(other: FunDBM[S, A])(implicit infField: InfField[A]): FunDBM[S, A]
  def decideState: DBMIxed[FunDBM, A]
  val innerMatrix: Option[FunMatrix[A]]
  def noOfVariables: Int
}

case class ClosedFunDBM[A](m: FunMatrix[A]) extends FunDBM[Closed, A] {
  require (m.dimension % 2 == 0)
  def noOfVariables : Int = (m.dimension / 2)
  override def liftFromInner(f: (FunMatrix[A]) => FunMatrix[A])
                            (implicit ifield: InfField[A]): FunDBM[Closed, A] =
    ClosedFunDBM(f(m))

  override def union(other: FunDBM[Closed, A])
                    (implicit infField: InfField[A]): FunDBM[Closed, A] = other match {
    case ClosedFunDBM(m2) => {
      val me: Matrix[FunMatrix] = FunMatrixMatrixInstance.funMatrixIsMatrix
      require(noOfVariables == other.noOfVariables)
      ClosedFunDBM(me.combine(infField.max)(m, m2))
    }
    case BottomFunDBM(n2) => { require(noOfVariables == other.noOfVariables) ; this }
  }

  val innerMatrix: Option[FunMatrix[A]] = Some(m)

  def decideState: DBMIxed[FunDBM, A] = CIxed(this)
}

case class NonClosedFunDBM[A](m: FunMatrix[A]) extends FunDBM[NonClosed, A] {
  require (m.dimension % 2 == 0)
  def noOfVariables : Int = (m.dimension / 2)
  override def liftFromInner(f: (FunMatrix[A]) => FunMatrix[A])
                            (implicit ifield: InfField[A]): FunDBM[NonClosed, A] =
    NonClosedFunDBM(f(m))
  def union(other: FunDBM[NonClosed, A])
           (implicit infField: InfField[A]): FunDBM[NonClosed, A] = other match {
    case NonClosedFunDBM(m2) => {
      val me: Matrix[FunMatrix] = FunMatrixMatrixInstance.funMatrixIsMatrix
      require(noOfVariables == other.noOfVariables)
      NonClosedFunDBM(me.combine(infField.max)(m, m2))
    }
  }

  val innerMatrix: Option[FunMatrix[A]] = Some(m)

  def decideState: DBMIxed[FunDBM, A] = NCIxed(this)
}

case class BottomFunDBM[A](noOfVariables: Int) extends FunDBM[Closed, A] {
  override def liftFromInner(f: (FunMatrix[A]) => FunMatrix[A])
                            (implicit ifield: InfField[A]): FunDBM[Closed, A] = {
    val bogus: FunMatrix[A] =
      FunMatrix((_, _) => ifield.infinity, noOfVariables * 2)
    BottomFunDBM(f(bogus).dimension / 2)
  }

  def union(other: FunDBM[Closed, A])(implicit infField: InfField[A]): FunDBM[Closed, A] =
    { require(noOfVariables == other.noOfVariables) ; other }

  val innerMatrix: Option[FunMatrix[A]] = None

  def decideState: DBMIxed[FunDBM, A] = {
    val dbm: FunDBM[Closed, A] = BottomFunDBM[A](noOfVariables)
    CIxed(dbm)
  }
}

object FunDBMInstance {
  val me: Matrix[FunMatrix] = FunMatrixMatrixInstance.funMatrixIsMatrix

  implicit val funDBM: DifferenceBoundMatrix[FunDBM] { type PosetConstraint[A] = InfField[A] } =
    new DifferenceBoundMatrix[FunDBM] {
    def update[S <: DBMState, A](f: (Int, Int) => A)
                                (m: FunDBM[S, A])
                                (implicit ifield: InfField[A]): ExistsM[A] =
      mkExFun(m.liftFromInner(me.update(f)))

    def incrementalClosure[S <: DBMState, A](v: VarIndex)
                             (dbm: FunDBM[S, A])
      (implicit evidence: InfField[A]): FunDBM[Closed, A] =
      dbm.innerMatrix match {
        case Some(m) =>
          BagnaraStrongClosure.incrementalClosure(v)(m) match {
            case Some(closed: FunMatrix[A]) => ClosedFunDBM(closed)
            case None => BottomFunDBM(dbm.noOfVariables)
          }
        case None => BottomFunDBM(dbm.noOfVariables)
      }

    def strongClosure[S <: DBMState, A](dbm: FunDBM[S, A])
                        (implicit evidence: InfField[A]): FunDBM[Closed, A] =
      dbm.innerMatrix match {
        case Some(m) =>
          BagnaraStrongClosure.strongClosure(m) match {
            case Some(closed: FunMatrix[A]) => ClosedFunDBM(closed)
            case None => BottomFunDBM(dbm.noOfVariables)
          }
        case None => BottomFunDBM(dbm.noOfVariables)
      }

    def forget[S <: DBMState, A](vi: VarIndex)(m: FunDBM[S, A])
                                (implicit ifield: InfField[A]): FunDBM[S, A] = {
      m.liftFromInner((inner) => {
        val f: (Int, Int) => A = (i, j) => {
          if (i != varPlus(vi) && i != varMinus(vi) && j != varPlus(vi) && j != varMinus(vi))
            me.get(i, j)(inner)
          else if (i == j && (i == varPlus(vi) || i == varMinus(vi)))
            ifield.zero
          else ifield.infinity
        }
        me.update(f)(inner)
      })
    }

    def nOfVars[S <: DBMState, A](m: FunDBM[S, A]): Int = m.noOfVariables

    def get[S <: DBMState, A](i: Int, j: Int)(m: FunDBM[S, A])
                             (implicit ifield: InfField[A]): Option[A] =
      m.innerMatrix.map((mmm) => me.get(i, j)(mmm))

    def mkExFun[S <: DBMState, A](funDBM: FunDBM[S, A]): ExistsM[A] =
      MkEx[S, ({ type T[S] = FunDBM[S, A]})#T](funDBM)

    def dbmIntersection[A, S <: DBMState, T <: DBMState]
      (m1: FunDBM[S, A], m2: FunDBM[T, A])
      (implicit ifield: InfField[A]): ExistsM[A] = {
      require(m1.noOfVariables == m2.noOfVariables)
      val o = for {
        mm1 <- m1.innerMatrix
        mm2 <- m2.innerMatrix
      } yield NonClosedFunDBM(me.combine(ifield.min)(mm1, mm2))
      o match {
        case Some(matrix) => mkExFun(matrix)
        case None => mkExFun(BottomFunDBM(m1.noOfVariables))
      }
    }

    def topDBM[A](nOfVars: Int)(implicit ifield: InfField[A]): FunDBM[Closed, A] =
      ClosedFunDBM(FunMatrix((i, j) => ifield.infinity, nOfVars * 2))

    def bottomDBM[A](nOfVars: Int)(implicit ifield: InfField[A]): FunDBM[Closed, A] = BottomFunDBM(nOfVars)
    def fromFun[A](d: Int, f: ((Int, Int) => A))(implicit ifield: InfField[A]): FunDBM[Closed, A] =
      strongClosure(NonClosedFunDBM(FunMatrix[A](f, d)))
    def flipVar[S <: DBMState, A](vi: VarIndex)(dbm: FunDBM[S, A])
                                 (implicit ifield: InfField[A]): FunDBM[S, A] = {
      dbm.liftFromInner((inner) => {
        val f: (Int, Int) => A = (i, j) => {
          if (i == varPlus(vi) || i == varMinus(vi)) {
            if (j == varPlus(vi) || j == varMinus(vi))
              me.get(signed(i), signed(j))(inner)
            else
              me.get(signed(i), j)(inner)
          } else {
            if (j == varPlus(vi) || j == varMinus(vi))
              me.get(i, signed(j))(inner)
            else
              me.get(i, j)(inner)
          }
        }
        me.update(f)(inner)
      })
    }

    def dbmUnion[S <: DBMState, A](m1: FunDBM[S, A], m2: FunDBM[S, A])
                                  (implicit ifield: InfField[A]): FunDBM[S, A] = m1.union(m2)

    def addScalarOnVar[S <: DBMState, A](vi: VarIndex, const: A)
                                        (fundbm: FunDBM[S, A])
                                        (implicit ifield: InfField[A]): FunDBM[S, A] = {
      fundbm.liftFromInner((dbm) => {
        val f: (Int, Int) => A = (i, j) => {
          val g1 = (i == varPlus(vi) && j != varPlus(vi) && j != varMinus(vi)) ||
                   (j == varMinus(vi) && i != varPlus(vi) && i != varMinus(vi))
          val g2 = (i != varPlus(vi) && i != varMinus(vi) && j == varPlus(vi)) ||
                   (j != varPlus(vi) && j != varMinus(vi) && i == varMinus(vi))
          val g3 = i == varPlus(vi) && j == varMinus(vi)
          val g4 = i == varMinus(vi) && j == varPlus(vi)
          if (g1) ifield.-(me.get(i, j)(dbm), const) else
          if (g2) ifield.+(me.get(i, j)(dbm), const) else
          if (g3) ifield.-(me.get(i, j)(dbm), ifield.double(const)) else
          if (g4) ifield.+(me.get(i, j)(dbm), ifield.double(const)) else
            me.get(i, j)(dbm)
        }
        me.update(f)(dbm)
      })
    }

    def isBottomDBM[A, S <: DBMState](m: FunDBM[S, A]): Boolean = m match {
      case BottomFunDBM(_) => true
      case _ => false
    }

    def widening[A, S <: DBMState, T <: DBMState]
      (dbm1: FunDBM[S, A], dbm2: FunDBM[T, A])
      (implicit ifield: InfField[A]): ExistsM[A] = {

      require(dbm1.noOfVariables == dbm2.noOfVariables)

      val m: Option[FunMatrix[A]] = for {
        m1 <- dbm1.innerMatrix
        m2 <- dbm2.innerMatrix
      } yield me.combine((mij: A, nij: A) =>
        ifield.compare(mij, nij) match {
          case GT => mij
          case _ => ifield.infinity
        }
      )(m1, m2)

      m match {
        case Some(matrix) => mkExFun(NonClosedFunDBM(matrix))
        case None => mkExFun(BottomFunDBM(dbm1.noOfVariables))
      }
    }

    def narrowing[A, S <: DBMState, T <: DBMState]
      (dbm1: FunDBM[S, A], dbm2: FunDBM[T, A])
      (implicit ifield: InfField[A]): ExistsM[A] = {

      require(dbm1.noOfVariables == dbm2.noOfVariables)

      val m = for {
        m1 <- dbm1.innerMatrix
        m2 <- dbm2.innerMatrix
      } yield me.combine((mij: A, nij: A) => if (mij == ifield.infinity) nij else mij)(m1, m2)
      m match {
        case Some(matrix) => mkExFun(NonClosedFunDBM(matrix))
        case None => mkExFun(BottomFunDBM(dbm1.noOfVariables))
      }
    }

    def decideState[S <: DBMState, A](dbm: FunDBM[S, A]): DBMIxed[FunDBM, A] =
      dbm.decideState

    def isTopDBM[A, S <: DBMState](dbm: FunDBM[S, A])
                                  (implicit ifield: InfField[A]): Boolean =
      dbm.innerMatrix match {
        case Some(m) =>
          me.toList(m)
            .foldRight(true)((el, b) => el == ifield.infinity && b)
        case None => false
      }

    def addVariable[S <: DBMState, A](dbm: FunDBM[S, A])
                                     (implicit ifield: InfField[A]): FunDBM[S, A] =
      dbm.liftFromInner((m) => {
        FunMatrix((i, j) => {
          if (i < dbm.noOfVariables * 2 && j < dbm.noOfVariables * 2)
            m(i, j)
          else ifield.infinity
        }, dbm.noOfVariables + 1)
      })

    // Proved with pen and paper that shuffling variables around preserves
    // closure state.
    def mapVariables[S <: DBMState, A](f: (VarIndex) => Option[VarIndex])
                                      (dbm: FunDBM[S, A])
                                      (implicit ifield: InfField[A]): FunDBM[S, A] =
      dbm.liftFromInner(VarMapping.mapVariablesAux(f, dbm.noOfVariables))

    def deleteVariable[S <: DBMState, A](v: VarIndex)(dbm: FunDBM[S, A])
                                        (implicit ifield: InfField[A]): FunDBM[S, A] =
      mapVariables((x : VarIndex) =>
        if (x.i < v.i) Some(x) else
        if (x.i == v.i) None else Some(VarIndex(x.i - 1)))(dbm)

    type PosetConstraint[A] = InfField[A]

    def compare[A](x: ExistsDBM[({ type Q[S] = FunDBM[S, A]})#Q],
                   y: ExistsDBM[({ type Q[S] = FunDBM[S, A]})#Q])
                  (implicit evidence: InfField[A]): Option[Ordering] = {
      (x.elem.innerMatrix, y.elem.innerMatrix) match {
        case (None, None) => Some(EQ)
        case (Some(_), None) => None
        case (None, Some(_)) => None
        case (Some(m1), Some(m2)) => {
          val l: List[Ordering] = me.combine(evidence.compare)(m1, m2).toList
          (l.forall(_ == EQ), l.forall(_ == LT), l.forall(_ == GT)) match {
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
object BagnaraStrongClosure {
  private val me: Matrix[FunMatrix] = FunMatrixMatrixInstance.funMatrixIsMatrix

  private def signed(i: Int) = if (i % 2 == 0) i + 1 else i - 1

  def nullCheck[A](m: FunMatrix[A])(implicit ifield: InfField[A]): Option[FunMatrix[A]] = {
    val negative: Boolean = (0 until m.dimension).exists((i) =>
      ifield.compare(me.get(i, i)(m), ifield.zero) == LT)
    if (negative) None else {
      val updater: (Int, Int) => A = (i, j) =>
        if (i == j) ifield.zero else me.get(i, j)(m)
      Some(me.update(updater)(m))
    }
  }

  def strengthen[A](dbm: FunMatrix[A])(implicit ifield: InfField[A]): FunMatrix[A] =
    (for { i <- 0 until dbm.dimension ;
           j <- 0 until dbm.dimension } yield (i, j))
      .foldLeft(dbm)((x, pair) => pair match {
        case (i, j) => {
          val newVal =
            ifield.min(
              me.get(i, j)(x),
              ifield.half(ifield.+(me.get(i, signed(i))(x), me.get(signed(j), j)(x))))
          me.update(i, j, newVal)(x)
        }
      })

  def strongClosure[A](dbm: FunMatrix[A])(implicit ifield: InfField[A]): Option[FunMatrix[A]] = {
    val closed: FunMatrix[A] = (for {
      k <- 0 until dbm.dimension
      i <- 0 until dbm.dimension
      j <- 0 until dbm.dimension
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

  def incrementalClosure[A](vi: VarIndex)(dbm: FunMatrix[A])
                           (implicit ifield: InfField[A]): Option[FunMatrix[A]] = {
    val p: ((Int, Int, Int)) => Boolean = {
      case (k, i, j) =>
        k == vi.i || k == signed(vi.i) ||
          i == vi.i || i == signed(vi.i) ||
          j == vi.i || j == signed(vi.i)
    }

    val iclosed: FunMatrix[A] =
      (for { k <- 0 until dbm.dimension ;
             i <- 0 until dbm.dimension ;
             j <- 0 until dbm.dimension } yield (k, i, j))
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

object VarMapping {

  def varMapImageSize(f: VarIndex => Option[VarIndex], nOfVars: Int): Int =
    (0 until nOfVars).map(VarIndex).map((v) => f(v) match {
      case Some(_) => 1
      case None => 0
    }).sum

  def mapVariablesAux[A](f: (VarIndex) => Option[VarIndex], nOfVars: Int)
                        (m: FunMatrix[A]): FunMatrix[A] = {

    def inverse(f: VarIndex => Option[VarIndex])(v: VarIndex): VarIndex = {
      val vars: List[Int] = (0 until nOfVars).toList
      val opts: List[Option[VarIndex]] = vars.map((n) => f(VarIndex(n)) match {
        case Some(vres) => ???
        case None => None
      })
      val opt: Option[VarIndex] = opts.flatten.headOption
      opt match {
        case Some(vres) => vres
        case None => throw new IllegalArgumentException("inverse out of bounds")
      }
    }

    def inverseIx(f: VarIndex => Option[VarIndex])(i: Int): Int =
      toIndexAndCoeff(i) match {
        case (v, Positive) => varPlus(inverse(f)(v))
        case (v, Negative) => varMinus(inverse(f)(v))
      }

    val newDim: Int = varMapImageSize(f, nOfVars)
    FunMatrix((i, j) => {
      m(inverseIx(f)(i), inverseIx(f)(j))
    }, newDim * 2)
  }
}
