package it.unich.jandom.domains.numerical.octagon.fast

import it.unich.jandom.domains.numerical.octagon._
import variables._
import VarIndexOps._
import CountOps._

case class HalfSubMatrix[A](mat: HalfMatrix[A], indices: Seq[VarIndex])

object HalfMatrixDenseSparseInstance {
  val sparseThreshold = 0.5

    val halfMatrixDenseSparseInstance = new DenseSparse[HalfMatrix] {

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
            DenseStrongClosure(m)

        def incrementalClosure[A](v: VarIndex)(m: HalfMatrix[A])
                                 (implicit e: InfField[A])
                                 : Option[HalfMatrix[A]] =
          if (
            variables.Fast.sparsityIndex(
              dimToVarCount(m.dimension),
              computeSparsity(m))
              >=  sparseThreshold)
            DenseIncrementalClosure(m, v)
          else
            SparseIncrementalClosure(m, v)

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

  val halfMatrixDecomposableInstance = new Decomposable[HalfMatrix, HalfSubMatrix] {

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

  def computeSparsity[A](m: HalfMatrix[A])(implicit ifield: InfField[A]): NNI =
    NNI(m.toSeq.count(v => ifield.compare(ifield.infinity, v) == EQ))
}

object DenseIncrementalClosure extends DenseClosureStrategy {
  def apply[A](m: HM[A], v: VarIndex)(implicit ifield: InfField[A]): Option[HM[A]] = {
    def br(v: VarIndex, k: VarIndex) = if (k.i < v.i) 2*k.i else 2*v.i
    // TODO: Adapt this into something `var`-less?
    var mutableM = m
    for (k <- allVars(e.nOfVars(m))) { // 0 to N
      for (i <- (2*v.i until 2*v.i + 2)) {
        val ik : A = e.get(i, 2*k.i)(mutableM)
        val ikk : A = e.get(i, 2*k.i + 1)(mutableM)
        // TODO We can merge these into a single loop, since e.get takes care of everything?
        for (j <- 0 until br(v, k)) {
          val kj : A = e.get(2*k.i, j)(mutableM)
          val kkj : A = e.get(2*k.i+1, j)(mutableM)
          mutableM = e.update(i, j,
            ifield.min(e.get(i,j)(m),
              ifield.min(ifield.+(ik, kj), ifield.+(ikk, kkj))
            )
          )(mutableM)
        }
        for (j <- br(v, k) until 2*v.i) {
          val kj : A = e.get(2*k.i, j)(mutableM)
          val kkj : A = e.get(2*k.i+1, j)(mutableM)
          mutableM = e.update(i, j,
            ifield.min(e.get(i,j)(m),
              ifield.min(ifield.+(ik, kj), ifield.+(ikk, kkj))
            )
          )(mutableM)
        }
      }

      for (j <- 2*v.i until 2*v.i + 2) {
        val kj : A = e.get(2*k.i, j)(mutableM)
        val kkj : A = e.get(2*k.i+1, j)(mutableM)
        for (i <-allIndices(varCountToDim(e.nOfVars(m))).drop(2*v.i)) { // 2v to n
          val ik : A = e.get(i, 2*k.i)(mutableM)
          val ikk : A = e.get(i, 2*k.i + 1)(mutableM)
          mutableM = e.update(i, j,
            ifield.min(e.get(i,j)(m),
              ifield.min(ifield.+(ik, kj), ifield.+(ikk, kkj))
            )
          )(mutableM)
        }
      }
    }
    val newM = loopBody(mutableM, v)
    nullCheck(strengtheningHalfScalar(newM))
  }
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

object SparseStrongClosure extends SparseClosureStrategy {
  def apply[A](m: HalfMatrix[A])
              (implicit ifield: InfField[A])
              : Option[HalfMatrix[A]] = {
    val newMat = allVars(e.nOfVars(m)).foldLeft(m)(loopBody)
    strengthening(newMat)
  }
}

object SparseIncrementalClosure extends SparseClosureStrategy {
  def apply[A](m: HalfMatrix[A], v: VarIndex)
    (implicit ifield: InfField[A])
      : Option[HalfMatrix[A]] = {
    def br(v: VarIndex, k: VarIndex) = if (k.i < v.i) 2*k.i else 2*v.i
    // TODO: Adapt this into something `var`-less?
    var mutableM = m
    for (k <- allVars(e.nOfVars(m))) { // 0 to N
      for (i <- (2*v.i until 2*v.i + 2)) {
        val ik : A = e.get(i, 2*k.i)(mutableM)
        if (ik != ifield.infinity) {
          for (j <- 0 until 2*v.i) {
            val kj : A = e.get(2*k.i, j)(mutableM)
            mutableM = e.update(i, j,
              ifield.min(e.get(i,j)(m),
                ifield.+(ik, kj))
            )(mutableM)
          }
        }
        val ikk : A = e.get(i, 2*k.i + 1)(mutableM)
        if (ikk != ifield.infinity) {
          for (j <- 0 until 2*v.i) {
            val kkj : A = e.get(2*k.i+1, j)(mutableM)
            mutableM = e.update(i, j,
              ifield.min(e.get(i,j)(m),
                ifield.+(ikk, kkj))
            )(mutableM)
          }
        }
      }

      for (j <- (2*v.i until 2*v.i + 2)) {
        val kj : A = e.get(2*k.i, j)(mutableM)
        if (kj != ifield.infinity) {
          for (i <- allIndices(varCountToDim(e.nOfVars(m))).drop(2*v.i)) { // 2v to n
            val kj : A = e.get(2*k.i, j)(mutableM)
            val ik : A = e.get(i, 2*k.i)(mutableM)
            mutableM = e.update(i, j,
              ifield.min(e.get(i,j)(m),
                ifield.+(ik, kj))
            )(mutableM)
          }
        }
        val kkj : A = e.get(2*k.i+1, j)(mutableM)
        if (kkj != ifield.infinity) {
          for (i <- allIndices(varCountToDim(e.nOfVars(m))).drop(j)) { // j to n
            val kkj : A = e.get(2*k.i+1, j)(mutableM)
            val ikk : A = e.get(i, 2*k.i + 1)(mutableM)
            mutableM = e.update(i, j,
              ifield.min(e.get(i,j)(m),
                ifield.+(ikk, kkj))
            )(mutableM)
          }
        }
      }
    }
    val newMat = loopBody(mutableM, v)
    strengthening(newMat)
  }
}

trait SparseClosureStrategy {

  type HM[A] = HalfMatrix[A]

  val e: DenseSparse[HalfMatrix] =
    HalfMatrixDenseSparseInstance.halfMatrixDenseSparseInstance

  def id[A](x: A) = x

  // Computes the location of finite values for 2k and 2k+1 -th row and column.
  protected def computeIndex[A]
    (m: HalfMatrix[A], k: Int)(implicit ifield: InfField[A])
    : (Seq[Int], Seq[Int], Seq[Int], Seq[Int]) = {

    val idxs = allIndices(varCountToDim(e.nOfVars(m)))

    val cp = idxs.drop(2 * k + 2).filter(i => m(i, 2*k) != ifield.infinity)
    val cm = idxs.drop(2 * k + 2).filter(i => m(i, 2*k + 1) != ifield.infinity)
    val rp = idxs.take(2 * k).filter(j => m(2*k, j) != ifield.infinity)
    val rm = idxs.take(2 * k).filter(j => m(2*k + 1, j) != ifield.infinity)

    (rp, rm, cp, cm)
  }

  // Uses the column indices to update the elements in j and j+1 -th column.
  protected def computeColumn[A]
    (m: HalfMatrix[A], j: Int, jj: Int, c: Seq[Int], cm: Seq[Int])
    (implicit ifield: InfField[A]) = {

    // update both m and c
    val (newM, newC) =
      if (m(j, jj) != ifield.infinity)
        cm.foldLeft((m, c))((pair, i) => {
          val (m, c) = pair
          val min = ifield.min(m(i, j), ifield.+(m(i, jj), m(jj, j)))
          val newMat = m.update(i, j, min)
          if (m(i, j) != ifield.infinity) (newMat, c) else (newMat, c :+ i)
        }) else (m, c)

    val t = (allIndices(m.dimension)).map(i => newM(i, j))
    (newM, newC, t)
  }

  protected def computeRow[A]
    (m: HalfMatrix[A], k: Int, kk: Int, r: Seq[Int], rm: Seq[Int])
    (implicit ifield: InfField[A]) = {

    if (m(k, kk) != ifield.infinity)
      rm.foldLeft((m, r))((pair, j) => {
        val (m, r) = pair
        val min = ifield.min(m(k, j), ifield.+(m(k, kk), m(kk, j)))
        val newMat = m.update(k, j, min)
        if (m(k, j) != ifield.infinity) (newMat, r) else (newMat, r :+ j)
      })
    else
      (m, r)
  }

  // Uses the indices obtained after 2k and 2k+1 -th row and column to compute
  // the remaining elements.
  protected def computeIteration[A]
    (m: HalfMatrix[A], k: Int, cp: Seq[Int], cm: Seq[Int],
     rp: Seq[Int], rm: Seq[Int], a: Seq[A], b: Seq[A])
    (implicit ifield: InfField[A]) = {

    def buildLoop
      (ikBldr: (HM[A], Int) => A, coll1: Seq[Int], coll2: Seq[Int],
       kjBldr1: (HM[A], Int) => A, kjBldr2: (HM[A], Int) => A,
       lMod1: Int => Int, rMod1: Int => Int, lMod2: Int => Int, rMod2: Int => Int)
      (m: HM[A], i1: Int): HM[A] = {

      val ik = ikBldr(m, i1)
      val m2  = coll1.foldLeft(m)((m, j1) => {
        val kj = kjBldr1(m, j1) ; val i = lMod1(i1) ; val j = rMod1(j1)
        m.update(i, j, ifield.min(m(i, j), ifield.+(ik, kj))) })
      val m3 = coll2.foldLeft(m2)((m, j1) => {
        val kj = kjBldr2(m, j1) ; val i = lMod2(i1) ; val j = rMod2(j1)
        m.update(i, j, ifield.min(m(i, j), ifield.+(ik, kj))) })
      m3
    }

    val mat1 = rp.foldLeft(m)(
      buildLoop((m, i1) => m(2 * k, i1), rm, cp, (m, j1) => m(2 * k + 1, j1),
        (m, j1) => b(j1), signed, id, signed, signed))
    val mat2 = rm.foldLeft(mat1)(
      buildLoop((m, i1) => m(2 * k + 1, i1), rp, cm, (m, j1) => m(2 * k, j1),
        (m, j1) => a(j1), signed, id, signed, signed))
    val mat3 = cm.foldLeft(mat2)(
      buildLoop((m, i1) => m(i1, 2 * k + 1), rm, cp, (m, j1) => m(2 * k + 1, j1),
        (m, j1) => b(j1), id, id, id, signed))
    val mat4 = cp.foldLeft(mat3)(
      buildLoop((m, i1) => m(i1, 2 * k), rp, cm, (m, j1) => m(2 * k, j1),
        (m, j1) => a(j1), id, id, id, signed))
    mat4
  }


  protected def strengthening[A](m: HalfMatrix[A])
                              (implicit ifield: InfField[A])
                              : Option[HalfMatrix[A]] = {
    val idxs = allIndices(varCountToDim(e.nOfVars(m)))
    val d = idxs.filter(i => m(signed(i), i) != ifield.infinity)

    val newMat: HalfMatrix[A] =
      d.foldLeft(m)((m, i) => {
        val ii = m(signed(i), i)
        d.foldLeft(m)((m, j) => {
          val jj = m(signed(j), j)
          val min = ifield.min(m(signed(i), j), ifield.half(ifield.+(ii,jj)))
          m.update(signed(i), j, min) }) })

    if (idxs.exists(i => ifield.compare(newMat(i, i), ifield.zero) == LT))
      None
    else
      Some(newMat)
  }

  def loopBody[A] (m: HalfMatrix[A], vi: VarIndex)(implicit ifield: InfField[A]) = {
    val (rp, rm, cp, cm) = computeIndex(m, vi.i)
    val (m1, cp1, a) = computeColumn(m, 2*vi.i, 2*vi.i + 1, cp, cm)
    val (m2, cm1, b) = computeColumn(m1, 2*vi.i + 1, 2*vi.i, cm, cp1)
    val (m3, rp1) = computeRow(m2, 2*vi.i, 2*vi.i + 1, rp, rm)
    val (m4, rm1) = computeRow(m3, 2*vi.i + 1, 2*vi.i, rm, rp1)
    computeIteration(m4, vi.i, rp1, rm1, cp1, cm1, a, b)
  }
}
