package it.unich.jandom.domains.numerical.octagon
// import scalaz.{Applicative, Monoid}

// A DBM is a matrix for which is makes sense to compute a strong closure.
trait DifferenceBoundMatrix[M[_]] extends Lattice1[M] {
  // strong closure and incremental closure are assumed to test for emptiness,
  // and return the bottom element in the positive case.
  def strongClosure[A](m: M[A])(implicit evidence: InfField[A]): M[A]
  def incrementalClosure[A](m: M[A])(implicit evidence: InfField[A]): M[A]
  def bottomDBM[A]: M[A]
  def topDBM[A]: M[A]
  def get[A](i: Int, j: Int)(m: M[A]): A
}

// // Provide some default implementations for DBMs over a field.
// trait FieldDBM[M[_]] extends DifferenceBoundMatrix[M] {

//   override type PosetConstraint[A] = InfField[A]
//   override def compare[A](x: M[A], y: M[A])(implicit evidence: PosetConstraint[A]): Option[Ordering] = {
//     val m: List[Option[Ordering]] = toList(combine(evidence.compare _)(x, y))
//     Applicative[Option].sequence(m) match {
//       case None => None
//       case Some(list) =>
//         (list.forall(_ == LT()), list.forall(_ == EQ()), list.forall(_ == GT())) match {
//           case (true, _, _) => Some(LT())
//           case (_, true, _) => Some(EQ())
//           case (_, _, true) => Some (GT())
//           case _ => None
//         }
//     }
//   }

//   override type LatticeConstraint[A] = InfField[A]
//   override def union[A](x: M[A], y: M[A])(implicit evidence: LatticeConstraint[A]): M[A] =
//     combine(evidence.max)(x, y)
//   override def intersection[A](x: M[A], y: M[A])(implicit evidence: InfField[A]): M[A] =
//     combine(evidence.min)(x, y)

// }
