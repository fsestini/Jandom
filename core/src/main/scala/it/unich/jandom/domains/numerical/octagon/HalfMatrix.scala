package it.unich.jandom.domains.numerical.octagon

/**
  * Created by fsestini on 7/11/17.
  *
  * Implementation of matrices that are represented by storing only the lower
  * triangular half of it, as explained in Singh et al.
  */
class HalfMatrix[A] {

}

object HalfMatrixIsRawDBM {
  implicit val halfMatrixIsRawDBM: RawDBM[HalfMatrix] = ???
}
