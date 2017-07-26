package it.unich.jandom.domains.numerical.octagon

import it.unich.jandom.domains.numerical.LinearForm
import it.unich.jandom.domains.numerical.octagon.VarIndexOps.OctaVarCoeff
import it.unich.jandom.utils.numberext.RationalExt
import it.unich.jandom.domains.numerical.octagon.VarIndexOps._
import spire.math.Rational

object OctagonalConstraint {

  sealed trait LFPosition
  case object ConstPos extends LFPosition
  case class VarPos(varIndex: VarIndex) extends LFPosition

  def buildLF(f: LFPosition => Option[Rational], dimension: Int): LinearForm = {
    def zeroIfNone(x : Option[Rational]): Rational =
      x match { case None => 0 ; case Some(y) => y }
    LinearForm(
      zeroIfNone(f(ConstPos)) ::
        List.range(0, dimension).map(VarIndex).map((vi) => zeroIfNone(f(vi))))
  }

  trait OctaConstraint[A]
  case class SingleConstraint[A](v: VarIndex, coeff: OctaVarCoeff,
                                 const: A, ifield: InfField[A])
    extends OctaConstraint[A]
  // The tacitly assume that v1 != v2
  case class DoubleConstraint[A](v1: VarIndex, coeff1: OctaVarCoeff,
                               v2: VarIndex, coeff2: OctaVarCoeff,
                               const: A, ifield: InfField[A])
    extends OctaConstraint[A]

  // TODO: Quite naive... We may want to avoid printing equivalent constraints,
  // TODO: like x - y <= c and - y + x <= c
  def prettyConstraint[A](p: A => String, vars: Seq[String])
                         (constr: OctaConstraint[A]): String =
    constr match {
      case SingleConstraint(v, coeff, const, ifield) => coeff match {
        case Positive => "v" + vars(v) + " <= " + p(const)
        case Negative => "-v" + vars(v) + " <= " + p(const)
      }
      case DoubleConstraint(v1, c1, v2, c2, c, ifield) => (c1, c2) match {
        case (Positive, Positive) => "v" + vars(v1) + " -" + vars(v2) + " <= " + p(c)
        case (Positive, Negative) => "v" + vars(v1) + " +" + vars(v2) + " <= " + p(c)
        case (Negative, Positive) => "-v" + vars(v1) + " -" + vars(v2) + " <= " + p(c)
        case (Negative, Negative) => "-v" + vars(v1) + " +" + vars(v2) + " <= " + p(c)
      }
    }

  def mapConstraint[A, B](f: A => B)
                         (c: OctaConstraint[A])
                         (implicit ifield: InfField[B]): OctaConstraint[B] =
    c match {
      case SingleConstraint(v, coeff, const, e) =>
        SingleConstraint(v, coeff, f(const), ifield)
      case DoubleConstraint(v1, c1, v2, c2, c, e) =>
        DoubleConstraint(v1, c1, v2, c2, f(c), ifield)
    }

  def octaConstrAt[S <: DBMState, A, M[_, _]](i: Int, j: Int, dbm: M[S, A])
                                    (implicit ifield: InfField[A],
                                     e: DifferenceBoundMatrix[M]): Option[OctaConstraint[A]] =
    (toIndexAndCoeff(i), toIndexAndCoeff(j)) match {
      case ((v1, c1), (v2, c2)) =>
        e.get(i, j)(dbm).map((c) => OctaConstraint(v1, c1, v2, c2, c, ifield))
    }

  def constrToLf(dimension: Int)
                (constr: OctaConstraint[RationalExt]): Option[LinearForm] = {
    constr match {
      case OctaConstraint(v1, c1, v2, c2, c, ifield) => {
        val rc1 = vcAsNumeral(c1)
        val rc2 = vcAsNumeral(c2)

        if (c < RationalExt.PositiveInfinity)
          Some(buildLF({
            case ConstPos => Some(- c.value)
            case VarPos(v) => if (v1 == v2) {
              val csum = rc1 + rc2
              if (csum == 0) None else
                Some(- (c.value / csum))
            } else {
              if (v == v1) Some(rc1)
              else if (v == v2) Some(rc2)
              else None
            }
          }, dimension))
        else None
      }
    }
  }

}