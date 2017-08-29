/**
  * Copyright 2017 Davide Dal Bianco, Filippo Sestini, Tobia Tesan
  *
  * This file is part of JANDOM: JVM-based Analyzer for Numerical DOMains
  * JANDOM is free software: you can redistribute it and/or modify
  * it under the terms of the GNU General Public License as published by
  * the Free Software Foundation, either version 3 of the License, or
  * (at your option) any later version.
  *
  * JANDOM is distributed in the hope that it will be useful,
  * but WITHOUT ANY WARRANTY; without even the implied warranty of a
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  * GNU General Public License for more details.
  *
  * You should have received a copy of the GNU General Public License
  * along with JANDOM.  If not, see <http://www.gnu.org/licenses/>.
  */

package it.unich.jandom.domains.numerical.octagon.fast
import it.unich.jandom.domains.numerical._
import it.unich.jandom.domains.numerical.octagon._
import it.unich.jandom.domains.numerical.octagon.variables.VarCount
import InfField._
import org.scalatest.FunSuite
import it.unich.jandom.utils.numberext.RationalExt

class FastOctagonSuite extends FunSuite {

  val mev = MEvidence(
    HalfMatrixDenseSparseInstance.halfMatrixDecomposableInstance,
    HalfMatrixDenseSparseInstance.halfMatrixDenseSparseInstance,
    HalfMatrixDenseSparseInstance.halfMatrixDecomposableInstance)

  val inf = Double.PositiveInfinity
  val e = CFDBMInstance.instance(mev)
  type FastDBM[S, A] = CFastDBM[HalfMatrix, HalfSubMatrix, S, A]
  val box =  BoxRationalDomain()
  val oct = new OctagonDomain[FastDBM, RationalExt, BoxRationalDomain](e, box)
  type DOM = OctagonDomain[FastDBM, RationalExt, BoxRationalDomain]

  test("Sanity check for pure") {
    assert(mev.dec.pure(VarCount(1), RationalExt.PositiveInfinity).vec.size == 4)
    assert(mev.dec.pure(VarCount(1), RationalExt.PositiveInfinity).vec.forall(_ == RationalExt.PositiveInfinity))
    assert(mev.dec.pure(VarCount(2), RationalExt.PositiveInfinity).vec.size == 12)
    assert(mev.dec.pure(VarCount(2), 0).vec.forall(_ == 0))
  }

  test("T {v0 <- 123}.toInterval == [123,123]") {
    val top = e.topDBM[RationalExt](VarCount(1))
    val t = AbstractOctagon[DOM, FastDBM, RationalExt, BoxRationalDomain](top, oct, box, e)
    assert(t.dbm == top)
    val a = t.linearAssignment(0, LinearForm.c(123))
    assert(a.toInterval.high.head == 123 & a.toInterval.low.head == 123)
  }
}
