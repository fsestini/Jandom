package it.unich.jandom.domains.numerical.octagon.fast

import it.unich.jandom.domains.numerical.octagon._
import VarIndexOps._
import CountOps._

case class HalfMatrixDenseSparseDBM[A](mat: HalfMatrix[A],
                                       indices: Seq[VarIndex],
                                       dimension: VarCount)

object HalfMatrixDenseSparseInstance {
    val instance = new DenseSparseDBM[HalfMatrixDenseSparseDBM] {

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

        def get[A](i: Int, j: Int)(m: HalfMatrixDenseSparseDBM[A])
                  (implicit e: InfField[A]): Option[A] = Some(m.mat(i, j))

        def varIndices[A](m: HalfMatrixDenseSparseDBM[A]): Seq[VarIndex] = m.indices

        def update[A](f: (Int, Int) => A)
                     (m: HalfMatrixDenseSparseDBM[A])
                     : HalfMatrixDenseSparseDBM[A] = {
          val elemIndices = varsToIndices(m.indices)
          val newMat = elemIndices.foldLeft(m.mat)({ case (mat, (i, j)) =>
              mat.update(i, j, f(i, j))
            })
          HalfMatrixDenseSparseDBM(newMat, m.indices, m.dimension)
        }

        def dbmUnion[A](m1: HalfMatrixDenseSparseDBM[A],
                        m2: HalfMatrixDenseSparseDBM[A])
                       (implicit e: InfField[A]): HalfMatrixDenseSparseDBM[A] = {
          val f = (i: Int, j: Int) => e.max(m1.mat(i, j), m2.mat(i, j))
          update(f)(m1)
        }

        def dbmIntersection[A](m1: HalfMatrixDenseSparseDBM[A],
                               m2: HalfMatrixDenseSparseDBM[A])
                              (implicit e: InfField[A])
                              : HalfMatrixDenseSparseDBM[A] = {
          val f = (i: Int, j: Int) => e.min(m1.mat(i, j), m2.mat(i, j))
          update(f)(m1)
        }

        def widening[A](m1: HalfMatrixDenseSparseDBM[A],
                        m2: HalfMatrixDenseSparseDBM[A])
                       (implicit e: InfField[A])
                       : HalfMatrixDenseSparseDBM[A] = {
          val f = (i: Int, j: Int) => {
            e.compare(m1.mat(i, j), m2.mat(i, j)) match {
              case GT => m1.mat(i, j)
              case _ => e.infinity
            }
          }
          update(f)(m1)
        }

        def narrowing[A](m1: HalfMatrixDenseSparseDBM[A],
                         m2: HalfMatrixDenseSparseDBM[A])
                       (implicit e: InfField[A])
                       : HalfMatrixDenseSparseDBM[A] = {
          val f = (i: Int, j: Int) =>
              if (m1.mat(i, j) == e.infinity) m2.mat(i, j) else m1.mat(i, j)
          update(f)(m1)
        }

        def strongClosure[A](m: HalfMatrixDenseSparseDBM[A])
                         (implicit e: InfField[A])
                         : Option[HalfMatrixDenseSparseDBM[A]] =
          if (FastDbmUtils.nuffSparse(m.dimension, computeSparsity(m)))
            SparseStrongClosure(m)
          else
            denseStrongClosure(m)

        def incrementalClosure[A](v: VarIndex)(m: HalfMatrixDenseSparseDBM[A])
                                 (implicit e: InfField[A])
                                 : Option[HalfMatrixDenseSparseDBM[A]] =
          if (FastDbmUtils.nuffSparse(m.dimension, computeSparsity(m)))
            denseIncrementalClosure(v)(m)
          else
            sparseIncrementalClosure(v)(m)

        def forget[A](v: VarIndex)(m: HalfMatrixDenseSparseDBM[A])
                     (implicit e: InfField[A]): HalfMatrixDenseSparseDBM[A] = {
          val f = (i: Int, j: Int) =>
            if (toIndexAndCoeff(i)._1 == v || toIndexAndCoeff(j)._1 == v)
              if (i == j) e.zero else e.infinity
            else
              m.mat(i, j)
          update(f)(m)
        }

        def flipVar[A](v: VarIndex)(m: HalfMatrixDenseSparseDBM[A])
                      : HalfMatrixDenseSparseDBM[A] = {
          val f = (i: Int, j: Int) =>
            if (i == varPlus(v) || i == varMinus(v))
              if (j == varPlus(v) || j == varMinus(v))
                m.mat(signed(i), signed(j))
              else
                m.mat(signed(i), j)
            else
              if (j == varPlus(v) || j == varMinus(v))
                m.mat(i, signed(j))
              else
                m.mat(i, j)
          update(f)(m)
        }

        def addScalarOnVar[A](v: VarIndex, c: A)(m: HalfMatrixDenseSparseDBM[A])
                             (implicit ifield: InfField[A])
                             : HalfMatrixDenseSparseDBM[A] = {
          val f = (i: Int, j: Int) => {
            val g1 = (i == varPlus(v) && j != varPlus(v) && j != varMinus(v)) ||
                     (j == varMinus(v) && i != varPlus(v) && i != varMinus(v))
            val g2 = (i != varPlus(v) && i != varMinus(v) && j == varPlus(v)) ||
                     (j != varPlus(v) && j != varMinus(v) && i == varMinus(v))
            val g3 = i == varPlus(v) && j == varMinus(v)
            val g4 = i == varMinus(v) && j == varPlus(v)
            if (g1) ifield.-(m.mat(i, j), c) else
            if (g2) ifield.+(m.mat(i, j), c) else
            if (g3) ifield.-(m.mat(i, j), ifield.double(c)) else
            if (g4) ifield.+(m.mat(i, j), ifield.double(c)) else
              m.mat(i, j)
          }
          update(f)(m)
        }


        def addVariable[A](m: HalfMatrixDenseSparseDBM[A])
                          (implicit ifield: InfField[A])
                          : HalfMatrixDenseSparseDBM[A] = {
          val nOfVars = addOne(m.dimension)
          val newVar = VarIndex(nOfVars.count -1)
          val newMat = new HalfMatrix(doubledVarCount(nOfVars), ifield.infinity)
                            .update(varPlus(newVar),
                                    varPlus(newVar),
                                    ifield.zero)
                            .update(varPlus(newVar),
                                    varPlus(newVar),
                                    ifield.zero)
          pour(m)(HalfMatrixDenseSparseDBM(newMat,
                                           m.indices :+ newVar,
                                           nOfVars))
        }

        def deleteVariable[A](m: HalfMatrixDenseSparseDBM[A])
                             (implicit ifield: InfField[A])
                             : HalfMatrixDenseSparseDBM[A] = {
          val nOfVars = subOne(m.dimension)
          val remVar = VarIndex(nOfVars.count)
          val newMat = new HalfMatrix(doubledVarCount(nOfVars), ifield.infinity)
          val newHMat = HalfMatrixDenseSparseDBM(newMat,
                                                 m.indices.filter(_ != remVar),
                                                 nOfVars)
          val f = (i: Int, j: Int) => m.mat(i, j)
          update(f)(newHMat)
        }

        def mapVariables[A](f: VarIndex => Option[VarIndex])
                           (m: HalfMatrixDenseSparseDBM[A])
                           (implicit ifield: InfField[A])
                           : HalfMatrixDenseSparseDBM[A] = {
          val newSize = VarCount(allVars(m.dimension).count(f(_).isDefined))
          val newIndices = m.indices.map(f(_)).collect({
              case Some(vi) => vi
            })
          val newMat = new HalfMatrix(doubledVarCount(newSize), ifield.infinity)
          // g is the inverse of f
          val g = allVars(m.dimension).map(vi => f(vi) -> vi).collect({
              case (Some(vi), vj) => vi -> vj
            }).toMap
          val updater = (i: Int, j: Int) => {
            val (vi, si) = toIndexAndCoeff(i)
            val (vj, sj) = toIndexAndCoeff(j)
            val ii = fromIndexAndCoeff(g(vi), si)
            val jj = fromIndexAndCoeff(g(vj), sj)
            m.mat(ii, jj)
          }
          HalfMatrixDenseSparseDBM(newMat.update(updater),
                                   newIndices,
                                   newSize)
        }

        def compare[A](m1: HalfMatrixDenseSparseDBM[A],
                       m2: HalfMatrixDenseSparseDBM[A])
                      (implicit ifield: InfField[A]): Option[Ordering] = {
          val elemIndices = varsToIndices(m1.indices)
          val ord = elemIndices.map({ case (i, j) =>
              ifield.compare(m1.mat(i, j), m2.mat(i, j))
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

        //////////////////////////////////////////////////////////////////////////////

        def extract[A](is: Seq[VarIndex])(m: HalfMatrixDenseSparseDBM[A])
                      : HalfMatrixDenseSparseDBM[A] =
          HalfMatrixDenseSparseDBM(m.mat, is, m.dimension)

        def pour[A](source: HalfMatrixDenseSparseDBM[A])
                   (dest: HalfMatrixDenseSparseDBM[A])
                   : HalfMatrixDenseSparseDBM[A] = {
          val sourceElemIndices = varsToIndices(source.indices)
          val f = (i: Int, j: Int) =>
            if (sourceElemIndices.contains((i, j)))
              source.mat(i, j)
            else
              dest.mat(i, j)
          update(f)(dest)
        }

        def nOfVars[A](m: HalfMatrixDenseSparseDBM[A]): VarCount =
          VarCount(m.indices.size)

        def pure[A](d: VarCount, x: A): HalfMatrixDenseSparseDBM[A] = {
          val mat = new HalfMatrix(doubledVarCount(d), x)
          val indices = allVars(d)
          HalfMatrixDenseSparseDBM(mat, indices, d)
        }

    }
}

object HalfMatrixDenseSparseDBM {
  def computeSparsity[A](m: HalfMatrixDenseSparseDBM[A])
                        (implicit ifield: InfField[A]): NNI =
    NNI(m.mat.toSeq.count(v => ifield.compare(ifield.infinity, v) == EQ))

  def denseStrongClosure[A](m: HalfMatrixDenseSparseDBM[A])
                        (implicit ifield: InfField[A]) = ???
  def sparseStrongClosure[A](m: HalfMatrixDenseSparseDBM[A])
                         (implicit ifield: InfField[A]) = ???

  def denseIncrementalClosure[A](vi: VarIndex)(m: HalfMatrixDenseSparseDBM[A])
                        (implicit ifield: InfField[A]) = ???
  def sparseIncrementalClosure[A](vi: VarIndex)(m: HalfMatrixDenseSparseDBM[A])
                         (implicit ifield: InfField[A]) = ???
}

object SparseStrongClosure {

  // returns (r, r', c, c')
  private def computeIndex[A](m: HalfMatrixDenseSparseDBM[A], k: Int)
                             (implicit ifield: InfField[A])
                             : (Seq[Int], Seq[Int], Seq[Int], Seq[Int]) = {
    val cIndices = m.indices.filter(_ > VarIndex(k))
              .flatMap(vi => Seq(varPlus(vi), varMinus(vi)))
    val rIndices = m.indices.filter(_ < VarIndex(k))
              .flatMap(vi => Seq(varPlus(vi), varMinus(vi)))

    val cp = cIndices.filter(i => m.mat(i, 2*k) == ifield.infinity)
    val cm = cIndices.filter(i => m.mat(signed(i), 2*k + 1) == ifield.infinity)
    val rp = rIndices.filter(j => m.mat(2*k, j) == ifield.infinity)
    val rm = rIndices.filter(j => m.mat(2*k + 1, signed(j)) == ifield.infinity)


    (rp, rm, cp, cm)
  }

  // returns (m, cp, t)
  private def computeColumn[A](m: HalfMatrixDenseSparseDBM[A],
                               k: Int,
                               kk: Int,
                               cp: Seq[Int],
                               cm: Seq[Int])
                               (implicit ifield: InfField[A]) = {
    // update both m and cm
    val newVals =
      if (m.mat(k, kk) != ifield.infinity) {
        cm.foldLeft((m, cp))((pair, i) => {
          val (m, cp) = pair
          val min = ifield.min(m.mat(i, k), ifield.+(m.mat(i, kk), m.mat(kk, k)))
          val newMat = HalfMatrixDenseSparseDBM(m.mat.update(i, k, min),
                                                m.indices,
                                                m.dimension)
          if (m.mat(i, k) != ifield.infinity) {
            (newMat, cp)
          } else {
            (newMat, cp :+ i)
          }
        })
      } else {
        (m, cp)
      }
    val t = for (i <- allIndices(doubledVarCount(m.dimension))) yield newVals._1.mat(i, k)
    (newVals._1, newVals._2, t)
  }

  // returns (m, rp)
  private def computeRow[A](m: HalfMatrixDenseSparseDBM[A],
                            k: Int,
                            kk: Int,
                            rp: Seq[Int],
                            rm: Seq[Int])
                            (implicit ifield: InfField[A]) = {
    if (m.mat(k, kk) != ifield.infinity) {
      rm.foldLeft((m, rp))((pair, j) => {
        val (m, rp) = pair
        val min = ifield.min(m.mat(k, j), ifield.+(m.mat(k, kk), m.mat(kk, j)))
        val newMat = HalfMatrixDenseSparseDBM(m.mat.update(k, j, min),
                                              m.indices,
                                              m.dimension)
        if (m.mat(k, j) != ifield.infinity) {
          (newMat, rp)
        } else {
          (newMat, rp :+ j)
        }
      })
    } else {
      (m, rp)
    }
  }

  private def computeIteration[A](m: HalfMatrixDenseSparseDBM[A],
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
        val ik = m.mat(2*k, i)
        val zero =
          rm.foldLeft(m)((m, j) => {
            val kj = m.mat(2*k + 1, j)
            val min = ifield.min(m.mat(signed(i), j), ifield.+(ik, kj))
            HalfMatrixDenseSparseDBM(m.mat.update(signed(i), j, min),
                                     m.indices,
                                     m.dimension)
          })
        cp.foldLeft(zero)((m, j) => {
            val kj = b(j)
            val min = ifield.min(m.mat(signed(i), signed(j)), ifield.+(ik, kj))
            HalfMatrixDenseSparseDBM(m.mat.update(signed(i), signed(j), min),
                                     m.indices,
                                     m.dimension)
          })
      })

    val mat2 =
      rm.foldLeft(mat1)((m, i) => {
        val ikk = m.mat(2*k + 1, i)
        val zero =
          rp.foldLeft(m)((m, j) => {
            val kkj = m.mat(2*k, j)
            val min = ifield.min(m.mat(signed(i), j), ifield.+(ikk, kkj))
            HalfMatrixDenseSparseDBM(m.mat.update(signed(i), j, min),
                                     m.indices,
                                     m.dimension)
          })
        cm.foldLeft(zero)((m, j) => {
            val kkj = a(j)
            val min = ifield.min(m.mat(signed(i), signed(j)), ifield.+(ikk, kkj))
            HalfMatrixDenseSparseDBM(m.mat.update(signed(i), signed(j), min),
                                     m.indices,
                                     m.dimension)
          })
      })

    val mat3 =
      cm.foldLeft(mat2)((m, i) => {
        val ik = m.mat(i, 2*k + 1)
        val zero =
          rm.foldLeft(m)((m, j) => {
            val kj = m.mat(2*k+1, j)
            val min = ifield.min(m.mat(i, j), ifield.+(ik, kj))
            HalfMatrixDenseSparseDBM(m.mat.update(i, j, min),
                                     m.indices,
                                     m.dimension)
          })
        cp.foldLeft(zero)((m, j) => {
            val kj = b(j)
            val min = ifield.min(m.mat(i, signed(j)), ifield.+(ik, kj))
            HalfMatrixDenseSparseDBM(m.mat.update(i, signed(j), min),
                                     m.indices,
                                     m.dimension)
          })
      })

    val mat4 =
      cp.foldLeft(mat3)((m, i) => {
        val ikk = m.mat(i, 2*k)
        val zero =
          rp.foldLeft(m)((m, j) => {
            val kkj = m.mat(2*k, j)
            val min = ifield.min(m.mat(i, j), ifield.+(ikk, kkj))
            HalfMatrixDenseSparseDBM(m.mat.update(i, j, min),
                                     m.indices,
                                     m.dimension)
          })
        cm.foldLeft(zero)((m, j) => {
            val kkj = a(j)
            val min = ifield.min(m.mat(i, signed(j)), ifield.+(ikk, kkj))
            HalfMatrixDenseSparseDBM(m.mat.update(i, signed(j), min),
                                     m.indices,
                                     m.dimension)
          })
      })

    mat4
  }


  private def strengthening[A](m: HalfMatrixDenseSparseDBM[A])
                              (implicit ifield: InfField[A])
                              : Option[HalfMatrixDenseSparseDBM[A]] = {
    val indices = m.indices.flatMap(vi => Seq(varPlus(vi), varMinus(vi)))
    val d = indices.filter(i => m.mat(signed(i), i) != ifield.infinity)
    val t = indices.map(i => m.mat(signed(i), i))

    val newMat: HalfMatrixDenseSparseDBM[A] =
      d.foldLeft(m)((m, i) => {
        val ii = t(i)
        d.foldLeft(m)((m, j) => {
          val jj = t(j)
          val min = ifield.min(m.mat(signed(i), j),
                               ifield.half(ifield.+(ii,jj)))

          HalfMatrixDenseSparseDBM(m.mat.update(signed(i), j, min),
                                   m.indices,
                                   m.dimension)
        })
      })

    if (indices.exists(i =>
            ifield.compare(newMat.mat(i, i), ifield.zero) == LT))
      None
    else
      Some(newMat)
  }

  def apply[A](m: HalfMatrixDenseSparseDBM[A])
              (implicit ifield: InfField[A])
              : Option[HalfMatrixDenseSparseDBM[A]] = {
    val newMat =
      m.indices.foldLeft(m)((m, vi) => {
        val (rp, rm, cp, cm) = computeIndex(m, vi.i)
        val (m1, cp1, a) = computeColumn(m, 2*vi.i, 2*vi.i + 1, cp, cm)
        val (m2, cm1, b) = computeColumn(m1, 2*vi.i + 1, 2*vi.i, cm, cp1)
        val (m3, rp1) = computeRow(m2, 2*vi.i, 2*vi.i + 1, rp, rm)
        val (m4, rm1) = computeRow(m3, 2*vi.i + 1, 2*vi.i, rm, rp1)
        computeIteration(m4, vi.i, rp1, rm1, cp1, cm1, a, b)
      })

    strengthening(m)
  }
}
