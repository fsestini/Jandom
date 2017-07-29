package it.unich.jandom.domains.numerical.octagon.fast

import it.unich.jandom.domains.numerical.octagon._
import VarIndexOps._

case class HalfMatrixDenseSparseDBM[A](mat: HalfMatrix[A],
                                       indeces: Seq[VarIndex],
                                       dimension: Int)

object HalfMatrixDenseSparseInstance {
    val instance = new DenseSparseDBM[HalfMatrixDenseSparseDBM] {

        private def varsToIndeces(indeces: Seq[VarIndex]): Seq[(Int, Int)] = {
          val varIndeces = for (vi <- indeces;
                                vj <- indeces;
                                if vi >= vj) yield (vi, vj)

          varIndeces.flatMap({ case (vi, vj) =>
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

        def varIndices[A](m: HalfMatrixDenseSparseDBM[A]): Seq[VarIndex] = m.indeces

        def update[A](f: (Int, Int) => A)
                     (m: HalfMatrixDenseSparseDBM[A])
                     : HalfMatrixDenseSparseDBM[A] = {
          val elemIndeces = varsToIndeces(m.indeces)
          val newMat = elemIndeces.foldLeft(m.mat)({ case (mat, (i, j)) =>
              mat.update(i, j, f(i, j))
            })
          HalfMatrixDenseSparseDBM(newMat, m.indeces, m.dimension)
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
                         : Option[HalfMatrixDenseSparseDBM[A]] = ???
        def incrementalClosure[A](v: VarIndex)(m: HalfMatrixDenseSparseDBM[A])
                                 (implicit e: InfField[A])
                                 : Option[HalfMatrixDenseSparseDBM[A]] = ???

        def forget[A](v: VarIndex)(m: HalfMatrixDenseSparseDBM[A])
                     : HalfMatrixDenseSparseDBM[A] = ???

        def flipVar[A](v: VarIndex)(m: HalfMatrixDenseSparseDBM[A])
                      : HalfMatrixDenseSparseDBM[A] = ???

        def addScalarOnVar[A](v: VarIndex, c: A)(m: HalfMatrixDenseSparseDBM[A])
                             (implicit ifield: InfField[A])
                             : HalfMatrixDenseSparseDBM[A] = ???


        def addVariable[A](m: HalfMatrixDenseSparseDBM[A])
                          (implicit ifield: InfField[A])
                          : HalfMatrixDenseSparseDBM[A] = ???
        def deleteVariable[A](m: HalfMatrixDenseSparseDBM[A])
                             (implicit ifield: InfField[A])
                             : HalfMatrixDenseSparseDBM[A] = ???
        def mapVariables[A](f: VarIndex => Option[VarIndex])
                           (m: HalfMatrixDenseSparseDBM[A])
                           (implicit ifield: InfField[A])
                           : HalfMatrixDenseSparseDBM[A] = ???
        def compare[A](m1: HalfMatrixDenseSparseDBM[A],
                       m2: HalfMatrixDenseSparseDBM[A])
                      (implicit ifield: InfField[A]): Option[Ordering] = ???

        //////////////////////////////////////////////////////////////////////////////

        def extract[A](is: Seq[VarIndex])(m: HalfMatrixDenseSparseDBM[A])
                      : HalfMatrixDenseSparseDBM[A] = ???
        def pour[A](source: HalfMatrixDenseSparseDBM[A])
                   (dest: HalfMatrixDenseSparseDBM[A])
                   : HalfMatrixDenseSparseDBM[A] = ???

        def nOfVars[A](m: HalfMatrixDenseSparseDBM[A]): Int = ???

        def pure[A](d: Int, x: A): HalfMatrixDenseSparseDBM[A] = ???
    }
}
