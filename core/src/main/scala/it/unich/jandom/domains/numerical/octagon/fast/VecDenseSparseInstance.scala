package it.unich.jandom.domains.numerical.octagon.fast

import it.unich.jandom.domains.numerical.octagon._
import variables._
import VarIndexOps._
import CountOps._
import scala.language.reflectiveCalls

case class VecSubMatrix[A](mat: VecMatrix[A], indices: Seq[VarIndex])

object VecMatrixDenseSparseInstance {

  val instance = new DenseSparse[VecMatrix] {

    def get[A](i: Int, j: Int)(m: VecMatrix[A]): A = m(i, j)

    def update[A](f: (Int, Int) => A)(m: VecMatrix[A]) = m.update(f)
    def update[A](i: Int, j: Int, x: A)(m: VecMatrix[A]) = m.update(i, j, x)

    def dbmUnion[A](m1: VecMatrix[A], m2: VecMatrix[A])
      (implicit ifield: InfField[A]) = {
      require(m1.dimension == m2.dimension, "dimension mismatch")
      m1.update((i: Int, j: Int) => ifield.max(m1(i, j), m2(i, j)))
    }

    def dbmIntersection[A](m1: VecMatrix[A], m2: VecMatrix[A])
      (implicit ifield: InfField[A]) = {
      require(m1.dimension == m2.dimension, "dimension mismatch")
      m1.update((i: Int, j: Int) => ifield.min(m1(i, j), m2(i, j)))
    }

    def widening[A](m1: VecMatrix[A], m2: VecMatrix[A])
      (implicit ifield: InfField[A]) = {
      require(m1.dimension == m2.dimension, "dimension mismatch")
      val f = (i: Int, j: Int) => {
        ifield.compare(m1(i, j), m2(i, j)) match {
          case GT => m1(i, j)
          case _ => ifield.infinity
        }
      }
      m1.update(f)
    }

    def narrowing[A](m1: VecMatrix[A], m2: VecMatrix[A])
      (implicit ifield: InfField[A]) = {
      require(m1.dimension == m2.dimension, "dimension mismatch")
      val f = (i: Int, j: Int) =>
      if (m1(i, j) == ifield.infinity) m2(i, j) else m1(i, j)
      m1.update(f)
    }

    def forget[A](v: VarIndex)(m: VecMatrix[A])(implicit ifield: InfField[A]) = {
      val f = (i: Int, j: Int) =>
        if (toIndexAndCoeff(i)._1 == v || toIndexAndCoeff(j)._1 == v)
          if (i == j) ifield.zero else ifield.infinity
        else m(i, j)
      m.update(f)
    }

    def flipVar[A](v: VarIndex)(m: VecMatrix[A]) = {
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
      m.update(f)
    }

    def addScalarOnVar[A](v: VarIndex, c: A)(m: VecMatrix[A])
      (implicit ifield: InfField[A]) = {
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
      m.update(f)
    }

    def compare[A](m1: VecMatrix[A], m2: VecMatrix[A])
      (implicit ifield: InfField[A]): Option[Ordering] = {
      require(m1.dimension == m2.dimension, "dimension mismatch")

      val ord = grid(m1.dimension).map(
        { case (i, j) => ifield.compare(m1(i, j), m2(i, j))})

      lazy val lt = ord.forall(v => v == EQ || v == LT)
      lazy val eq = ord.forall(v => v == EQ)
      lazy val gt = ord.forall(v => v == EQ || v == GT)

      if (lt) Some(LT) else if (eq) Some(EQ) else if (gt) Some(GT) else None
    }


    def strongClosure[A](m: VecMatrix[A])
      (implicit e: InfField[A]): Option[VecMatrix[A]] =
      (new BagnaraStrongClosure[VecMatrix]()(
        VecMatrixMatrixInstance.vecMatrixIsMatrix)).strongClosure(m)

    def incrementalClosure[A](v: VarIndex)(m: VecMatrix[A])
      (implicit e: InfField[A]): Option[VecMatrix[A]] =
      (new BagnaraStrongClosure[VecMatrix]()(
        VecMatrixMatrixInstance.vecMatrixIsMatrix)).incrementalClosure(v)(m)

    def nOfVars[A](m: VecMatrix[A]): VarCount = dimToVarCount(m.dimension)

  }

}

object VecMatrixDecomposableInstance {

  val instance = new Decomposable[VecMatrix, VecSubMatrix] {

    val ds = VecMatrixDenseSparseInstance.instance

    def addVariable[A](m: VecMatrix[A])(implicit ifield: InfField[A]) = {
      VecMatrix(
        varCountToDim(addOne(ds.nOfVars(m))),
        (i: Int, j: Int) => {
          if (!isComprised(ds.nOfVars(m))(toIndexAndCoeff(i)._1)
            || !isComprised(ds.nOfVars(m))(toIndexAndCoeff(j)._1))
            ifield.infinity else m(i, j)
        })
    }

    def deleteVariable[A](m: VecMatrix[A])(implicit ifield: InfField[A]) =
      VecMatrix(varCountToDim(subOne(ds.nOfVars(m))), m(_, _))


    def mapVariables[A](f: VarIndex => Option[VarIndex])(m: VecMatrix[A])
      (implicit ifield: InfField[A]): VecMatrix[A] = {

      val newSize = VarCount(allVars(ds.nOfVars(m)).count(f(_).isDefined))
      val newIndices =
        allVars(ds.nOfVars(m)).map(v => f(v)).collect({ case Some(vi) => vi })
      val newMat = VecMatrix(varCountToDim(newSize), ifield.infinity)

      // g is the inverse of f
      val g = allVars(ds.nOfVars(m))
        .map(vi => f(vi) -> vi)
        .collect({ case (Some(vi), vj) => vi -> vj }).toMap
      val updater = (i: Int, j: Int) => {
        val (vi, si) = toIndexAndCoeff(i)
        val (vj, sj) = toIndexAndCoeff(j)
        val ii = fromIndexAndCoeff(g(vi), si)
        val jj = fromIndexAndCoeff(g(vj), sj)
        m(ii, jj)
      }
      newMat.update(updater)
    }

    def extract[A](is: Seq[VarIndex])(m: VecMatrix[A]): VecSubMatrix[A] = {
      val subDim = is.length * 2
      val f: (Int, Int) => A = (i, j) => {
        val (vi, signi) = toIndexAndCoeff(i)
        val originalVi = is(vi.i)
        val (vj, signj) = toIndexAndCoeff(j)
        val originalVj = is(vj.i)
        m(fromIndexAndCoeff(originalVi, signi),
          fromIndexAndCoeff(originalVj, signj))
      }
      VecSubMatrix(VecMatrix(varCountToDim(VarCount(is.length)), f), is)
    }

    // TODO: not extremely efficient as it is. consider improvements
    def pour[A](source: VecSubMatrix[A])(dest: VecMatrix[A]): VecMatrix[A] = {
      def findCorresponding(
        ixs: Seq[VarIndex], v: VarIndex, sign: OctaVarCoeff): Option[Int] =
        ixs.zipWithIndex
          .find(vyvx => vyvx._1 == v)
          .map(p => VarIndex(p._2))
          .map(v => fromIndexAndCoeff(v, sign))

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

    def pure[A](d: VarCount, x: A) = VecMatrix(varCountToDim(d), x)
    def varIndices[A](m: VecSubMatrix[A]): Seq[VarIndex] = m.indices
    def compare[A](m1: VecSubMatrix[A], m2: VecSubMatrix[A])
      (implicit ifield: InfField[A]) = ds.compare(m1.mat, m2.mat)
    def widening[A](m1: VecSubMatrix[A], m2: VecSubMatrix[A])
                   (implicit e: InfField[A]): VecSubMatrix[A] =
      VecSubMatrix(ds.widening(m1.mat, m2.mat), m1.indices)
    def update[A](f: (Int, Int) => A)(m: VecSubMatrix[A]): VecSubMatrix[A] =
      VecSubMatrix(ds.update(f)(m.mat), m.indices)

    def update[A](i: Int, j: Int, x: A)(m: VecSubMatrix[A]): VecSubMatrix[A] =
      VecSubMatrix(ds.update(i, j, x)(m.mat), m.indices)

    def strongClosure[A](m: VecSubMatrix[A])(implicit e: InfField[A]) =
      ds.strongClosure(m.mat).map(clo => VecSubMatrix(clo, m.indices))

    def get[A](i: Int, j: Int)(m: VecSubMatrix[A]): A = ds.get(i, j)(m.mat)

    def dbmIntersection[A](m1: VecSubMatrix[A], m2: VecSubMatrix[A])
      (implicit e: InfField[A]): VecSubMatrix[A] =
      VecSubMatrix(ds.dbmIntersection(m1.mat, m2.mat), m1.indices)

    def dbmUnion[A](m1: VecSubMatrix[A], m2: VecSubMatrix[A])
      (implicit e: InfField[A]): VecSubMatrix[A] =
      VecSubMatrix(ds.dbmUnion(m1.mat, m2.mat), m1.indices)

    def incrementalClosure[A](v: VarIndex)(m: VecSubMatrix[A])
      (implicit e: InfField[A]): Option[VecSubMatrix[A]] =
      ds.incrementalClosure(v)(m.mat).map(clo => VecSubMatrix(clo, m.indices))

    def forget[A](v: VarIndex)(m: VecSubMatrix[A])(implicit e: InfField[A]): VecSubMatrix[A] =
      VecSubMatrix(ds.forget(v)(m.mat), m.indices)

    def flipVar[A](v: VarIndex)(m: VecSubMatrix[A]): VecSubMatrix[A] =
      VecSubMatrix(ds.flipVar(v)(m.mat), m.indices)

    def narrowing[A](m1: VecSubMatrix[A], m2: VecSubMatrix[A])
      (implicit e: InfField[A]): VecSubMatrix[A] =
      VecSubMatrix(ds.narrowing(m1.mat, m2.mat), m1.indices)

    def addScalarOnVar[A](v: VarIndex, c: A)(m: VecSubMatrix[A])
      (implicit ifield: InfField[A]): VecSubMatrix[A] =
      VecSubMatrix(ds.addScalarOnVar(v, c)(m.mat), m.indices)

  }

}
