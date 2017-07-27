package it.unich.jandom.domains.numerical.octagon

import it.unich.jandom.domains.numerical.{DenseLinearForm, LinearForm}
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
    DenseLinearForm(
      zeroIfNone(f(ConstPos)) ::
        List.range(0, dimension).map(VarIndex).map((vi) => zeroIfNone(f(VarPos(vi)))))
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
        case Positive => "v" + vars(v.i) + " <= " + p(const)
        case Negative => "-v" + vars(v.i) + " <= " + p(const)
      }
      case DoubleConstraint(v1, c1, v2, c2, c, ifield) => (c1, c2) match {
        case (Positive, Positive) => "v" + vars(v1.i) + " -" + vars(v2.i) + " <= " + p(c)
        case (Positive, Negative) => "v" + vars(v1.i) + " +" + vars(v2.i) + " <= " + p(c)
        case (Negative, Positive) => "-v" + vars(v1.i) + " -" + vars(v2.i) + " <= " + p(c)
        case (Negative, Negative) => "-v" + vars(v1.i) + " +" + vars(v2.i) + " <= " + p(c)
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

  def octaConstrAt[S <: DBMState, A, M[_, _]]
    (i: Int, j: Int, dbm: M[S, A])
    (implicit ifield: InfField[A], e: DifferenceBoundMatrix[M]): Option[OctaConstraint[A]] =

    (toIndexAndCoeff(i), toIndexAndCoeff(j)) match {
      case ((v1, c1), (v2, c2)) =>
        e.get(i, j)(dbm).flatMap((c) => {
          val rc1: A = vcAsNumeral[A](c1)
          val rc2: A = vcAsNumeral[A](c2)

          if (c != RationalExt.PositiveInfinity) {
            if (v1 == v2) {
              val csum = ifield.-(rc1, rc2)

              ifield.compare(csum, ifield.zero) match {
                case EQ => None
                case GT => Some(SingleConstraint(v1, Positive, ifield./(c, csum), ifield))
                case LT => Some(SingleConstraint(v1, Negative, ifield./(c, ifield.inverse(csum)), ifield))
              }

            } else Some(DoubleConstraint(v1, c1, v2, c2, c, ifield))
          } else None
        })
    }

  def constrToLf(dimension: Int)
                (constr: OctaConstraint[RationalExt]): Option[LinearForm] =

    constr match {
      case SingleConstraint(v, coeff, const, ifield) =>
        if (const < RationalExt.PositiveInfinity) {
          Some(buildLF({
            case ConstPos => Some(- const.value)
            case VarPos(vv) => if (v == vv) Some(vcAsNumeral[Rational](coeff)) else None
          }, dimension))
        } else None
      case DoubleConstraint(v1, c1, v2, c2, c, ifield) =>
        if (c < RationalExt.PositiveInfinity) {
          Some(buildLF({
            case ConstPos => Some(- c.value)
            case VarPos(v) =>
              if (v == v1) Some(vcAsNumeral[Rational](c1))
              else if (v == v2) Some(- vcAsNumeral[Rational](c2))
              else None
          }, dimension))
        } else None
    }

}
