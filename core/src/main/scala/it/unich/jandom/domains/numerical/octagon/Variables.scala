package it.unich.jandom.domains.numerical.octagon

// Distinguish integers used as variable indices
case class VarIndex(i: Int) extends Ordered[VarIndex] {
  def compare(that: VarIndex) = this.i compare that.i
}
