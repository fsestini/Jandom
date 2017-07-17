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

package it.unich.jandom.domains.numerical

import org.scalatest.FunSuite

class OctagonDomainSuite extends FunSuite {
  lazy val dom = new OctagonDomain()
  import octagon._

  sealed abstract class DummyDbm[+A,B](){}

  case class RegularDummyDbm[+A, B]
    (val arr : Option[Seq[Seq[B]]],
      val state : A,
      val nOfVars : Int)
      extends DummyDbm[A,B]
  case class TopDummyDbm[B] (
    val nOfVars : Int) extends  DummyDbm[Closed ,B]
  case class BottomDummyDbm[B] (
    val nOfVars : Int) extends  DummyDbm[Closed ,B]

  val e : DifferenceBoundMatrix[DummyDbm] = new DifferenceBoundMatrix[DummyDbm] {
    def get[B](i: Int, j: Int)(m: DummyDbm[DBMState, B]): Option[B] = m match {
      case m : RegularDummyDbm[_,_] => Some(m.arr.get(i)(j))
      case _ => ???
    }

    def update[B](f: (Int, Int) => B)(m: DummyDbm[DBMState, B]): DummyDbm[DBMState, B] = m match {
      case m : RegularDummyDbm[DBMState,B] => {
        m.copy(arr =
          Some(
            (0 to m.arr.get.size-1).map(
              (i) => (0 to (m.arr.get)(i).size-1).map(
                (j) => f(i,j)
              )
            )
          )
        )
      }
      case _ => ???
    }

    def strongClosure[B](m: DummyDbm[DBMState, B])
      (implicit evidence: InfField[B]): DummyDbm[Closed, B] = m match {
      case m : RegularDummyDbm[_,_] => m.copy(state = Closed())
      // Unsound! For testing only
      case m : TopDummyDbm[B] => m
      case m : BottomDummyDbm[B] => m
    }

    def incrementalClosure[B](v: VarIndex)(m: DummyDbm[DBMState, B])
      (implicit evidence: InfField[B]): DummyDbm[Closed, B] =  ???

    def bottomDBM[B](nOfVars : Int) : DummyDbm[Closed, B] =
      new BottomDummyDbm[B](nOfVars)

    def isBottomDBM[B](m: DummyDbm[DBMState, B]): Boolean = m match {
      case m : BottomDummyDbm[B] => true
      case _ => false
    }

    def topDBM[B](nOfVars : Int): DummyDbm[Closed, B] =
      new TopDummyDbm[B](nOfVars)

    def dbmUnion[S <: DBMState, B](m1: DummyDbm[S, B], m2: DummyDbm[S, B])
      (implicit ifield: InfField[B]): DummyDbm[S, B] = ???

    def dbmIntersection[B](m1: DummyDbm[DBMState, B], m2: DummyDbm[DBMState, B])
      (implicit ifield: InfField[B]): DummyDbm[DBMState, B] = ???

    def flipVar[S <: DBMState, B](v: VarIndex)
      (m: DummyDbm[S, B])(implicit ifield: InfField[B]): DummyDbm[S, B] = ???

    def addScalarOnVar[S <: DBMState, B](v: VarIndex, c: B)
      (m: DummyDbm[S, B])(implicit ifield: InfField[B]): DummyDbm[S, B] = ???

    def forget[S <: DBMState, B](v: VarIndex)(m: DummyDbm[S, B]): DummyDbm[S, B]  = ???

    def nOfVars[B](m: DummyDbm[DBMState, B]): Int = m match {
      case m : RegularDummyDbm[_,_] => m.nOfVars
      case m : TopDummyDbm[_] => m.nOfVars
      case m : BottomDummyDbm[_] => m.nOfVars
    }
  }

  test("Test singleton interval") {
    val arr = Some(Seq(Seq(0,-246.0),Seq(246.0,0)))
    val foo = new RegularDummyDbm[Closed, Double](arr, Closed(), 1)
    val ao = AbstractOctagon[DummyDbm](foo, e)
    assert(ao.toInterval.low(0) == ao.toInterval.high(0))
    assert(ao.toInterval.low(0) == 123)
  }

  test("Test join") {
    val arr = Some(Seq(Seq(0,-220.0),Seq(220.0,0)))
    val arr2 = Some(Seq(Seq(0,-200.0),Seq(100.0,0)))
    val foo = new RegularDummyDbm[Closed, Double](arr, Closed(), 1)
    val foo2 = new RegularDummyDbm[Closed, Double](arr2, Closed(), 1)
    val ao = AbstractOctagon[DummyDbm](foo, e)
    val ao2 = AbstractOctagon[DummyDbm](foo2, e)
    assert(ao.toInterval.low(0) == 110)
    assert(ao.toInterval.high(0) == 110)
    assert(ao2.toInterval.low(0) == 100)
    assert(ao2.toInterval.high(0) == 100)
    val join = ao.join(ao2)
    assert(join.toInterval.low(0) == 100)
    assert(join.toInterval.high(0) == 110)
  }
}
