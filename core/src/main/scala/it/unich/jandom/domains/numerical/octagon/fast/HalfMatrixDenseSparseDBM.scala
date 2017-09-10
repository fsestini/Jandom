package it.unich.jandom.domains.numerical.octagon.fast

import it.unich.jandom.domains.numerical.octagon._
import variables._
import VarIndexOps._
import CountOps._

case class HalfSubMatrix[A](mat: HalfMatrix[A], indices: Seq[VarIndex])

object HalfMatrixDenseSparseInstance {
  val sparseThreshold = 0.5

    val halfMatrixDenseSparseInstance = new DenseSparse[HalfMatrix] {

        import HalfMatrixDenseSparseDBM._

        private def varsToIndices(indices: Seq[VarIndex]): Seq[(Int, Int)] = {
          val varIndices = for (vi <- indices;
                                vj <- indices;
                                if vi >= vj) yield (vi, vj)

          varIndices.flatMap({ case (vi, vj) =>
            Seq(
              (varPlus(vi), varPlus(vj)),
              (varPlus(vi), varMinus(vj)),
              (varMinus(vi), varPlus(vj)),
              (varMinus(vi), varMinus(vj))
            )
          })
        }

        def get[A](i: Int, j: Int)(m: HalfMatrix[A]): A = m(i, j)

        def update[A](f: (Int, Int) => A)(m: HalfMatrix[A]): HalfMatrix[A] = m.update(f)

        def dbmUnion[A](m1: HalfMatrix[A], m2: HalfMatrix[A])
                       (implicit e: InfField[A]): HalfMatrix[A] = {
          require(m1.dimension == m2.dimension)
          val f = (i: Int, j: Int) => e.max(m1(i, j), m2(i, j))
          HalfMatrix(f, nOfVars(m1))
        }

        def dbmIntersection[A](m1: HalfMatrix[A], m2: HalfMatrix[A])
                              (implicit e: InfField[A]): HalfMatrix[A] = {
          require(m1.dimension == m2.dimension)
          val f = (i: Int, j: Int) => e.min(m1(i, j), m2(i, j))
          HalfMatrix(f, nOfVars(m1))
        }

        def widening[A](m1: HalfMatrix[A], m2: HalfMatrix[A])
                       (implicit e: InfField[A]): HalfMatrix[A] = {
          val f = (i: Int, j: Int) => {
            e.compare(m1(i, j), m2(i, j)) match {
              case GT => m1(i, j)
              case _ => e.infinity
            }
          }
          update(f)(m1)
        }

        def narrowing[A](m1: HalfMatrix[A], m2: HalfMatrix[A])
                        (implicit e: InfField[A]): HalfMatrix[A] = {
          val f = (i: Int, j: Int) =>
              if (m1(i, j) == e.infinity) m2(i, j) else m1(i, j)
          update(f)(m1)
        }

        def strongClosure[A](m: HalfMatrix[A])
                            (implicit e: InfField[A]): Option[HalfMatrix[A]] =
          if (
            variables.Fast.sparsityIndex(
              dimToVarCount(m.dimension),
              computeSparsity(m))
              >=  sparseThreshold)
            SparseStrongClosure(m)
          else
            denseStrongClosure(m)

        def incrementalClosure[A](v: VarIndex)(m: HalfMatrix[A])
                                 (implicit e: InfField[A])
                                 : Option[HalfMatrix[A]] =
          if (
            variables.Fast.sparsityIndex(
              dimToVarCount(m.dimension),
              computeSparsity(m))
              >=  sparseThreshold)
            denseIncrementalClosure(v)(m)
          else
            sparseIncrementalClosure(v)(m)

        def forget[A](v: VarIndex)(m: HalfMatrix[A])
                     (implicit e: InfField[A]): HalfMatrix[A] = {
          val f = (i: Int, j: Int) =>
            if (toIndexAndCoeff(i)._1 == v || toIndexAndCoeff(j)._1 == v)
              if (i == j) e.zero else e.infinity
            else
              m(i, j)
          update(f)(m)
        }

        def flipVar[A](v: VarIndex)(m: HalfMatrix[A])
                      : HalfMatrix[A] = {
          val f = (i: Int, j: Int) =>
            if (i == varPlus(v) || i == varMinus(v))
              if (j == varPlus(v) || j == varMinus(v))
                m(signed(i), signed(j))
              else
                m(signed(i), j)
            else
              if (j == varPlus(v) || j == varMinus(v))
                m(i, signed(j))
              else
                m(i, j)
          update(f)(m)
        }

        def addScalarOnVar[A](v: VarIndex, c: A)(m: HalfMatrix[A])
                             (implicit ifield: InfField[A])
                             : HalfMatrix[A] = {
          val f = (i: Int, j: Int) => {
            val g1 = (i == varPlus(v) && j != varPlus(v) && j != varMinus(v)) ||
                     (j == varMinus(v) && i != varPlus(v) && i != varMinus(v))
            val g2 = (i != varPlus(v) && i != varMinus(v) && j == varPlus(v)) ||
                     (j != varPlus(v) && j != varMinus(v) && i == varMinus(v))
            val g3 = i == varPlus(v) && j == varMinus(v)
            val g4 = i == varMinus(v) && j == varPlus(v)
            if (g1) ifield.-(m(i, j), c) else
            if (g2) ifield.+(m(i, j), c) else
            if (g3) ifield.-(m(i, j), ifield.double(c)) else
            if (g4) ifield.+(m(i, j), ifield.double(c)) else
              m(i, j)
          }
          update(f)(m)
        }


        def compare[A](m1: HalfMatrix[A], m2: HalfMatrix[A])
                      (implicit ifield: InfField[A]): Option[Ordering] = {
          require(m1.dimension == m2.dimension)
          val ord = m1.lowerIndices.map({ case (i, j) =>
              ifield.compare(m1(i, j), m2(i, j))
            })
          lazy val lt = ord.forall(v => v == EQ || v == LT)
          lazy val eq = ord.forall(v => v == EQ)
          lazy val gt = ord.forall(v => v == EQ || v == GT)
          if (lt)
            Some(LT)
          else if (eq)
            Some(EQ)
          else if (gt)
            Some(GT)
          else
            None
        }

      def update[A](i: Int, j: Int, x: A)(m: HalfMatrix[A]): HalfMatrix[A] = m.update(i, j, x)

      def nOfVars[A](m: HalfMatrix[A]): VarCount = dimToVarCount(m.dimension)
    }

  val halfMatrixDecomposableInstance: Decomposable[HalfMatrix, HalfSubMatrix] = new Decomposable[HalfMatrix, HalfSubMatrix] {

    val ds = halfMatrixDenseSparseInstance

    def addVariable[A](m: HalfMatrix[A])
                      (implicit ifield: InfField[A]): HalfMatrix[A] = {
      HalfMatrix((i, j) => {
        if (!isComprised(ds.nOfVars(m))(toIndexAndCoeff(i)._1)
          || !isComprised(ds.nOfVars(m))(toIndexAndCoeff(j)._1))
          ifield.infinity else m(i, j)
      }, addOne(ds.nOfVars(m)))
    }

    def deleteVariable[A](m: HalfMatrix[A])
                         (implicit ifield: InfField[A]): HalfMatrix[A] =
      HalfMatrix(m(_, _), subOne(ds.nOfVars(m)))

    def mapVariables[A](f: VarIndex => Option[VarIndex])(m: HalfMatrix[A])
                       (implicit ifield: InfField[A]): HalfMatrix[A] = {
      val newSize = VarCount(allVars(ds.nOfVars(m)).count(f(_).isDefined))
      val newIndices = allVars(ds.nOfVars(m)).map(v => f(v)).collect({
        case Some(vi) => vi
      })
      val newMat = new HalfMatrix(varCountToDim(newSize), ifield.infinity)
      // g is the inverse of f
      val g = allVars(ds.nOfVars(m)).map(vi => f(vi) -> vi).collect({
        case (Some(vi), vj) => vi -> vj
      }).toMap
      val updater = (i: Int, j: Int) => {
        val (vi, si) = toIndexAndCoeff(i)
        val (vj, sj) = toIndexAndCoeff(j)
        val ii = fromIndexAndCoeff(g(vi), si)
        val jj = fromIndexAndCoeff(g(vj), sj)
        m(ii, jj)
      }
      newMat.update(updater)
    }

    def extract[A](is: Seq[VarIndex])(m: HalfMatrix[A]): HalfSubMatrix[A] = {
      val subDim = is.length * 2
      val f: (Int, Int) => A = (i, j) => {
        val (vi, signi) = toIndexAndCoeff(i)
        val originalVi = is(vi.i)
        val (vj, signj) = toIndexAndCoeff(j)
        val originalVj = is(vj.i)
        m(fromIndexAndCoeff(originalVi, signi),
          fromIndexAndCoeff(originalVj, signj))
      }
      HalfSubMatrix(HalfMatrix(f, VarCount(is.length)), is)
    }

    // TODO: not extremely efficient as it is. consider improvements
    def pour[A](source: HalfSubMatrix[A])(dest: HalfMatrix[A]): HalfMatrix[A] = {
      def findCorresponding(ixs: Seq[VarIndex], v: VarIndex, sign: OctaVarCoeff): Option[Int] = {
        ixs.zipWithIndex
          .find(vyvx => vyvx._1 == v)
          .map(p => VarIndex(p._2))
          .map(v => fromIndexAndCoeff(v, sign))
      }
      val f: (Int, Int) => A = (i, j) => {
        val (vi, signi) = toIndexAndCoeff(i)
        val (vj, signj) = toIndexAndCoeff(j)
        (findCorresponding(source.indices, vi, signi),
          findCorresponding(source.indices, vj, signj)) match {
          case (Some(ii), Some(jj)) => ds.get(ii, jj)(source.mat)
          case _ => ds.get(i, j)(dest)
        }
      }
      ds.update(f)(dest)
    }

    def pure[A](d: VarCount, x: A): HalfMatrix[A] =
      new HalfMatrix[A](varCountToDim(d), x)

    def varIndices[A](m: HalfSubMatrix[A]): Seq[VarIndex] = m.indices
    def compare[A](m1: HalfSubMatrix[A], m2: HalfSubMatrix[A])
                  (implicit ifield: InfField[A]): Option[Ordering] = ds.compare(m1.mat, m2.mat)

    def widening[A](m1: HalfSubMatrix[A], m2: HalfSubMatrix[A])
                   (implicit e: InfField[A]): HalfSubMatrix[A] =
      HalfSubMatrix(ds.widening(m1.mat, m2.mat), m1.indices)

    def update[A](f: (Int, Int) => A)(m: HalfSubMatrix[A]): HalfSubMatrix[A] =
      HalfSubMatrix(ds.update(f)(m.mat), m.indices)

    def update[A](i: Int, j: Int, x: A)(m: HalfSubMatrix[A]): HalfSubMatrix[A] =
      HalfSubMatrix(ds.update(i, j, x)(m.mat), m.indices)

    def strongClosure[A](m: HalfSubMatrix[A])(implicit e: InfField[A]): Option[HalfSubMatrix[A]] =
      ds.strongClosure(m.mat).map(clo => HalfSubMatrix(clo, m.indices))

    def get[A](i: Int, j: Int)(m: HalfSubMatrix[A]): A = ds.get(i, j)(m.mat)

    def dbmIntersection[A](m1: HalfSubMatrix[A], m2: HalfSubMatrix[A])(implicit e: InfField[A]): HalfSubMatrix[A] =
      HalfSubMatrix(ds.dbmIntersection(m1.mat, m2.mat), m1.indices)

    def dbmUnion[A](m1: HalfSubMatrix[A], m2: HalfSubMatrix[A])(implicit e: InfField[A]): HalfSubMatrix[A] =
      HalfSubMatrix(ds.dbmUnion(m1.mat, m2.mat), m1.indices)

    def incrementalClosure[A](v: VarIndex)(m: HalfSubMatrix[A])(implicit e: InfField[A]): Option[HalfSubMatrix[A]] =
      ds.incrementalClosure(v)(m.mat).map(clo => HalfSubMatrix(clo, m.indices))

    def forget[A](v: VarIndex)(m: HalfSubMatrix[A])(implicit e: InfField[A]): HalfSubMatrix[A] =
      HalfSubMatrix(ds.forget(v)(m.mat), m.indices)

    def flipVar[A](v: VarIndex)(m: HalfSubMatrix[A]): HalfSubMatrix[A] =
      HalfSubMatrix(ds.flipVar(v)(m.mat), m.indices)

    def narrowing[A](m1: HalfSubMatrix[A], m2: HalfSubMatrix[A])(implicit e: InfField[A]): HalfSubMatrix[A] =
      HalfSubMatrix(ds.narrowing(m1.mat, m2.mat), m1.indices)

    def addScalarOnVar[A](v: VarIndex, c: A)(m: HalfSubMatrix[A])(implicit ifield: InfField[A]): HalfSubMatrix[A] =
      HalfSubMatrix(ds.addScalarOnVar(v, c)(m.mat), m.indices)
  }

}

object HalfMatrixDenseSparseDBM {
  def computeSparsity[A](m: HalfMatrix[A])
                        (implicit ifield: InfField[A]): NNI =
    NNI(m.toSeq.count(v => ifield.compare(ifield.infinity, v) == EQ))

  def denseStrongClosure[A](m: HalfMatrix[A])(implicit ifield: InfField[A]) =
    DenseStrongClosure(m)

  def sparseStrongClosure[A](m: HalfMatrix[A])(implicit ifield: InfField[A]) =
    SparseStrongClosure.apply(m)

  def denseIncrementalClosure[A](vi: VarIndex)(m: HalfMatrix[A])(implicit ifield: InfField[A]) =
    DenseStrongClosure(m)
  def sparseIncrementalClosure[A](vi: VarIndex)(m: HalfMatrix[A])(implicit ifield: InfField[A]) =
    SparseStrongClosure.apply(m)
}

object DenseIncrementalClosure {
  // TODO stub
  def apply[A](m: HalfMatrix[A])(implicit ifield: InfField[A]) =
    DenseStrongClosure(m)(ifield)
}

object SparseIncrementalClosure {
  // TODO stub
  def apply[A](m: HalfMatrix[A])(implicit ifield: InfField[A]) =
    SparseStrongClosure.apply(m)(ifield)
}

object DenseStrongClosure extends DenseClosureStrategy {
  def apply[A](m: HM[A])(implicit ifield: InfField[A]): Option[HM[A]] = {
    val newM = allVars(e.nOfVars(m)).foldLeft(m)(loopBody)
    nullCheck(strengtheningHalfScalar(newM))
  }
}

// Strong closure for dense half matrices.
// Taken from Singh, Fast Algorithms for Octagon Abstract Domain
trait DenseClosureStrategy {

  type HM[A] = HalfMatrix[A]

  protected val dec: Decomposable[HalfMatrix, HalfSubMatrix] =
    HalfMatrixDenseSparseInstance.halfMatrixDecomposableInstance
  protected val e: DenseSparse[HalfMatrix] =
    HalfMatrixDenseSparseInstance.halfMatrixDenseSparseInstance

  protected def computeColHalfScalar[A]
    (c: Int, d: Int)(m: HM[A])(implicit ifield: InfField[A]): HM[A] = {

    val s: Int = if (c % 2 == 1) c + 1 else c + 2
    val kj: A = e.get(d, c)(m)

    allIndices(varCountToDim(e.nOfVars(m))).drop(s). // s until n
      foldLeft(m)((mm, i) => {
      val ik: A = e.get(i, d)(mm)
      e.update(i, c, ifield.min(e.get(i,c)(mm), ifield.+(ik, kj)))(mm)
    })
  }

  protected def computeRowHalfScalar[A]
    (r: Int, s: Int)(m: HM[A])(implicit ifield: InfField[A]): HM[A] = {

    val ee: Int = if (r % 2 == 1) r - 1 else r
    val ik: A = e.get(r, s)(m)

    (0 until ee).foldLeft(m)((mm, j) => {
      val kj: A = e.get(s, j)(mm)
      e.update(r, j, ifield.min(e.get(r, j)(mm), ifield.+(ik, kj)))(mm)
    })
  }

  protected def computeIterationHalfScalar[A]
    (k: Int)(m: HM[A])(implicit ifield: InfField[A]): HM[A] = {

    def loop(m: HM[A], i: Int)(implicit ifield: InfField[A]): HM[A] = {
      val i2: Int = if (i % 2 == 0) i + 2 else i + 1
      val br = if (i2 < 2 * k) i2 else 2 * k
      val ik: A = e.get(i, 2 * k)(m)
      val ikk = e.get(i, 2 * k + 1)(m)
      val mbr = (0 until br).foldLeft(m)((m, j) => {
        val kj = e.get(2 * k, j)(m)
        val kkj = e.get(2 * k + 1, j)(m)
        e.update(i, j,
          ifield.min(e.get(i, j)(m),
            ifield.min(ifield.+(ik, kj), ifield.+(ikk, kkj))))(m)
      })
      (2 * k + 2 until i2).foldLeft(mbr)((m, j) => {
        val kj: A = e.get(j, 2 * k)(m) // a_j
        val kkj: A = e.get(j, 2 * k + 1)(m) // b_j
        e.update(i, j,
          ifield.min(e.get(i, j)(m),
            ifield.min(ifield.+(ik, kj), ifield.+(ikk, kkj))))(m)
      })
    }
    val indices = allIndices(varCountToDim(e.nOfVars(m)))
    val first = indices.take(2 * k). // 0 until 2k
      foldLeft(m)(loop)
    val second = indices.drop(2 * k + 2). // 2k+2 until n
      foldLeft(first)(loop)
    second
  }

  def strengtheningHalfScalar[A](m : HM[A])(implicit ifield: InfField[A]): HM[A] = {
    allIndices(varCountToDim(e.nOfVars(m))) // 0 until n
      .foldLeft(m)((m, i) => {
      val i2 : Int = if (i % 2 == 0) i + 2 else i + 1
      val ii: A = e.get(i, signed(i))(m)
      (0 until i2).foldLeft(m)((m, j) => {
        val jj = e.get(signed(j), j)(m)
        e.update(i, j,
          ifield.min(e.get(i, j)(m),
            ifield.half(ifield.+(ii, jj))))(m)
      })
    })
  }

  protected def nullCheck[A](m : HM[A])(implicit ifield: InfField[A]): Option[HM[A]] = {
    val bottom: Boolean = allIndices(varCountToDim(e.nOfVars(m))) // 0 until n
      .exists(
      i => ifield.compare(e.get(i, i)(m), ifield.zero) == LT)
    if (bottom) None else Some(m)
  }

  def loopBody[A](m: HM[A], k: VarIndex)(implicit ifield: InfField[A]) = {
    val f1 = computeColHalfScalar[A](2 * k.i, 2 * k.i + 1) _
    val f2 = computeColHalfScalar[A](2 * k.i + 1, 2 * k.i) _
    val f3 = computeRowHalfScalar[A](2 * k.i, 2 * k.i + 1) _
    val f4 = computeRowHalfScalar[A](2 * k.i + 1, 2 * k.i) _
    val f5 = computeIterationHalfScalar[A](k.i) _
    (f1 andThen f2 andThen f3 andThen f4 andThen f5)(m)
  }
}

object SparseStrongClosure {

  val e: DenseSparse[HalfMatrix] = HalfMatrixDenseSparseInstance.halfMatrixDenseSparseInstance

  def indices[A](m:HalfMatrix[A]): Seq[VarIndex] = allVars(e.nOfVars(m))

  // returns (r, r', c, c')
  private def computeIndex[A](m: HalfMatrix[A], k: Int)
                             (implicit ifield: InfField[A])
                             : (Seq[Int], Seq[Int], Seq[Int], Seq[Int]) = {
    val cIndices = indices(m).filter(_ > VarIndex(k))
              .flatMap(vi => Seq(varPlus(vi), varMinus(vi)))
    val rIndices = indices(m).filter(_ < VarIndex(k))
              .flatMap(vi => Seq(varPlus(vi), varMinus(vi)))

    val cp = cIndices.filter(i => m(i, 2*k) == ifield.infinity)
    val cm = cIndices.filter(i => m(signed(i), 2*k + 1) == ifield.infinity)
    val rp = rIndices.filter(j => m(2*k, j) == ifield.infinity)
    val rm = rIndices.filter(j => m(2*k + 1, signed(j)) == ifield.infinity)

    (rp, rm, cp, cm)
  }

  // returns (m, cp, t)
  private def computeColumn[A](m: HalfMatrix[A],
                               k: Int,
                               kk: Int,
                               cp: Seq[Int],
                               cm: Seq[Int])
                               (implicit ifield: InfField[A]) = {
    // update both m and cm
    val newVals =
      if (m(k, kk) != ifield.infinity) {
        cm.foldLeft((m, cp))((pair, i) => {
          val (m, cp) = pair
          val min = ifield.min(m(i, k), ifield.+(m(i, kk), m(kk, k)))
          val newMat = m.update(i, k, min)
          if (m(i, k) != ifield.infinity) {
            (newMat, cp)
          } else {
            (newMat, cp :+ i)
          }
        })
      } else {
        (m, cp)
      }
    val t = for (i <- allIndices(m.dimension)) yield newVals._1(i, k)
    (newVals._1, newVals._2, t)
  }

  // returns (m, rp)
  private def computeRow[A](m: HalfMatrix[A],
                            k: Int,
                            kk: Int,
                            rp: Seq[Int],
                            rm: Seq[Int])
                            (implicit ifield: InfField[A]) = {
    if (m(k, kk) != ifield.infinity) {
      rm.foldLeft((m, rp))((pair, j) => {
        val (m, rp) = pair
        val min = ifield.min(m(k, j), ifield.+(m(k, kk), m(kk, j)))
        val newMat = m.update(k, j, min)
        if (m(k, j) != ifield.infinity) {
          (newMat, rp)
        } else {
          (newMat, rp :+ j)
        }
      })
    } else {
      (m, rp)
    }
  }

  private def computeIteration[A](m: HalfMatrix[A],
                                  k: Int,
                                  cp: Seq[Int],
                                  cm: Seq[Int],
                                  rp: Seq[Int],
                                  rm: Seq[Int],
                                  a: Seq[A],
                                  b: Seq[A])
                                 (implicit ifield: InfField[A]) = {
    val mat1 =
      rp.foldLeft(m)((m, i) => {
        val ik = m(2*k, i)
        val zero =
          rm.foldLeft(m)((m, j) => {
            val kj = m(2*k + 1, j)
            val min = ifield.min(m(signed(i), j), ifield.+(ik, kj))
            m.update(signed(i), j, min)
          })
        cp.foldLeft(zero)((m, j) => {
            val kj = b(j)
            val min = ifield.min(m(signed(i), signed(j)), ifield.+(ik, kj))
            m.update(signed(i), signed(j), min)
          })
      })

    val mat2 =
      rm.foldLeft(mat1)((m, i) => {
        val ikk = m(2*k + 1, i)
        val zero =
          rp.foldLeft(m)((m, j) => {
            val kkj = m(2*k, j)
            val min = ifield.min(m(signed(i), j), ifield.+(ikk, kkj))
            m.update(signed(i), j, min)
          })
        cm.foldLeft(zero)((m, j) => {
            val kkj = a(j)
            val min = ifield.min(m(signed(i), signed(j)), ifield.+(ikk, kkj))
            m.update(signed(i), signed(j), min)
          })
      })

    val mat3 =
      cm.foldLeft(mat2)((m, i) => {
        val ik = m(i, 2*k + 1)
        val zero =
          rm.foldLeft(m)((m, j) => {
            val kj = m(2*k+1, j)
            val min = ifield.min(m(i, j), ifield.+(ik, kj))
            m.update(i, j, min)
          })
        cp.foldLeft(zero)((m, j) => {
            val kj = b(j)
            val min = ifield.min(m(i, signed(j)), ifield.+(ik, kj))
            m.update(i, signed(j), min)
          })
      })

    val mat4 =
      cp.foldLeft(mat3)((m, i) => {
        val ikk = m(i, 2*k)
        val zero =
          rp.foldLeft(m)((m, j) => {
            val kkj = m(2*k, j)
            val min = ifield.min(m(i, j), ifield.+(ikk, kkj))
            m.update(i, j, min)
          })
        cm.foldLeft(zero)((m, j) => {
            val kkj = a(j)
            val min = ifield.min(m(i, signed(j)), ifield.+(ikk, kkj))
            m.update(i, signed(j), min)
          })
      })

    mat4
  }


  private def strengthening[A](m: HalfMatrix[A])
                              (implicit ifield: InfField[A])
                              : Option[HalfMatrix[A]] = {
    val indicess = indices(m).flatMap(vi => Seq(varPlus(vi), varMinus(vi)))
    val d = indicess.filter(i => m(signed(i), i) != ifield.infinity)
    val t = indicess.map(i => m(signed(i), i))

    val newMat: HalfMatrix[A] =
      d.foldLeft(m)((m, i) => {
        val ii = t(i)
        d.foldLeft(m)((m, j) => {
          val jj = t(j)
          val min = ifield.min(m(signed(i), j),
                               ifield.half(ifield.+(ii,jj)))

          m.update(signed(i), j, min)
        })
      })

    if (indicess.exists(i =>
            ifield.compare(newMat(i, i), ifield.zero) == LT))
      None
    else
      Some(newMat)
  }

  def apply[A](m: HalfMatrix[A])
              (implicit ifield: InfField[A])
              : Option[HalfMatrix[A]] = {
    val newMat =
      indices(m).foldLeft(m)((m, vi) => {
        val (rp, rm, cp, cm) = computeIndex(m, vi.i)
        val (m1, cp1, a) = computeColumn(m, 2*vi.i, 2*vi.i + 1, cp, cm)
        val (m2, cm1, b) = computeColumn(m1, 2*vi.i + 1, 2*vi.i, cm, cp1)
        val (m3, rp1) = computeRow(m2, 2*vi.i, 2*vi.i + 1, rp, rm)
        val (m4, rm1) = computeRow(m3, 2*vi.i + 1, 2*vi.i, rm, rp1)
        computeIteration(m4, vi.i, rp1, rm1, cp1, cm1, a, b)
      })

    strengthening(newMat)
  }
}
