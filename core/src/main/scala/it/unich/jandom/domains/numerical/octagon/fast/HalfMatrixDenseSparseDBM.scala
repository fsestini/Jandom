package it.unich.jandom.domains.numerical.octagon.fast

import it.unich.jandom.domains.numerical.octagon._


case class HalfMatrixDenseSparseDBM[A](mat: HalfMatrix[A],
                                       indeces: Seq[VarIndex],
                                       dimension: Int)

object HalfMatrixDenseSparseInstance {
    val instance = new DenseSparseDBM[HalfMatrixDenseSparseDBM] {

        def get[A](i: Int, j: Int)(m: HalfMatrixDenseSparseDBM[A])
                  (implicit e: InfField[A]): Option[A] = Some(m.mat(i, j))

        def varIndices[A](m: HalfMatrixDenseSparseDBM[A]): Seq[VarIndex] = ???

        def update[A](f: (Int, Int) => A)
                     (m: HalfMatrixDenseSparseDBM[A])
                     : HalfMatrixDenseSparseDBM[A] = ???

        def dbmUnion[A](m1: HalfMatrixDenseSparseDBM[A],
                        m2: HalfMatrixDenseSparseDBM[A])
                       (implicit e: InfField[A]): HalfMatrixDenseSparseDBM[A] = ???
        def dbmIntersection[A](m1: HalfMatrixDenseSparseDBM[A],
                               m2: HalfMatrixDenseSparseDBM[A])
                              (implicit e: InfField[A])
                              : HalfMatrixDenseSparseDBM[A] = ???

        def widening[A](m1: HalfMatrixDenseSparseDBM[A],
                        m2: HalfMatrixDenseSparseDBM[A])
                       (implicit e: InfField[A])
                       : HalfMatrixDenseSparseDBM[A] = ???
        def narrowing[A](m1: HalfMatrixDenseSparseDBM[A],
                         m2: HalfMatrixDenseSparseDBM[A])
                       (implicit e: InfField[A])
                       : HalfMatrixDenseSparseDBM[A] = ???

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
