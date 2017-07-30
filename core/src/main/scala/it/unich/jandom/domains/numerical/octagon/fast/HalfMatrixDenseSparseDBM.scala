package it.unich.jandom.domains.numerical.octagon.fast

import it.unich.jandom.domains.numerical.octagon._
import VarIndexOps._

case class HalfMatrixDenseSparseDBM[A](mat: HalfMatrix[A],
                                       indices: Seq[VarIndex],
                                       dimension: Int)

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
            sparseStrongClosure(m)
          else
            denseStrongClosure(m)

        def incrementalClosure[A](v: VarIndex)(m: HalfMatrixDenseSparseDBM[A])
                                 (implicit e: InfField[A])
                                 : Option[HalfMatrixDenseSparseDBM[A]] = ???

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
          val nOfVars = m.mat.dimension + 1
          val newVar = VarIndex(nOfVars-1)
          val newMat = new HalfMatrix(nOfVars, ifield.infinity)
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
          val nOfVars = m.mat.dimension - 1
          val remVar = VarIndex(nOfVars)
          val newMat = new HalfMatrix(nOfVars, ifield.infinity)
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
          val allVarIndices = for (i <- 0 until m.dimension) yield VarIndex(i)
          val newSize = allVarIndices.count(f(_).isDefined)
          val newIndices = m.indices.map(f(_)).collect({
              case Some(vi) => vi
            })
          val newMat = new HalfMatrix(newSize, ifield.infinity)
          // g is the inverse of f
          val g = allVarIndices.map(vi => f(vi) -> vi).collect({
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

        def nOfVars[A](m: HalfMatrixDenseSparseDBM[A]): Int = m.indices.size

        def pure[A](d: Int, x: A): HalfMatrixDenseSparseDBM[A] = {
          val mat = new HalfMatrix(d, x)
          val indices = 0 until d map (VarIndex(_))
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
}



