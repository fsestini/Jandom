/**
 * Copyright 2013 Gianluca Amato
 * 
 * This file is part of JANDOM: JVM-based Analyzer for Numerical DOMains
 * JANDOM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * JANDOM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty ofa
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with JANDOM.  If not, see <http://www.gnu.org/licenses/>.
 */

package it.unich.sci.jandom
package gui

import domains.{ NumericalDomain, NumericalProperty }
import scala.swing._
import java.awt.GridBagConstraints
import scala.swing.ListView.Renderer
import javax.swing.SpinnerNumberModel
import javax.swing.JSpinner
import it.unich.sci.jandom.targets.{Parameters,Target}
import it.unich.sci.jandom.parameters.Parameter
import it.unich.sci.jandom.parameters.ParameterEnumeration
import it.unich.sci.jandom.parameters.ParameterValue
import it.unich.sci.jandom.parameters.WideningScope
import it.unich.sci.jandom.parameters.NarrowingStrategy
import it.unich.sci.jandom.widenings.DelayedWidening
import it.unich.sci.jandom.widenings.DefaultWidening
import it.unich.sci.jandom.narrowings.DelayedNarrowing
import it.unich.sci.jandom.narrowings.NoNarrowing
import it.unich.sci.jandom.widenings.DelayedWideningFactory
import it.unich.sci.jandom.narrowings.DelayedNarrowingFactory

class ParametersPane extends GridBagPanel {
  border = Swing.EmptyBorder(5, 5, 5, 5)
  val domainComboBox = addParameterEnumeration(0, NumericalDomain)
  val wideningComboBox = addParameterEnumeration(1, WideningScope)
  val narrowingComboBox = addParameterEnumeration(2, NarrowingStrategy)
  val delayModel = new SpinnerNumberModel(0, 0, Double.PositiveInfinity, 1)
  val delay = Component.wrap(new JSpinner(delayModel))
  val debug = new CheckBox("Debug")
  
  layout(new Label("Widening Delay:")) = new Constraints(0, 3, 1, 1, 0.0, 0.0, GridBagConstraints.EAST,
      GridBagConstraints.NONE, new Insets(0,0,5,5), 0, 0)
  layout(delay) = new Constraints(1, 3, 1, 1, 0.0, 0.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets(0,0,5,0), 0, 0)
  layout(debug) = new Constraints(0, 4, 2, 1, 0.0, 0.0, GridBagConstraints.BASELINE,
      GridBagConstraints.HORIZONTAL, new Insets(0,0,0,0), 0, 0)
  layout(Swing.VGlue) = new Constraints(0, 5, 2, 1, 0.0, 1.0, GridBagConstraints.BASELINE,
      GridBagConstraints.NONE, new Insets(0,0,0,0), 0, 0)
  
   object ParameterRenderer extends Renderer[ParameterValue] {
    val r = implicitly[Renderer[String]]
    def componentFor(list: ListView[_], isSelected: Boolean,
      focused: Boolean, a: ParameterValue, index: Int): Component =
      {
        val c = r.componentFor(list, isSelected, focused, a.name, index)
        c.tooltip = a.description
        return c
      }
  }
  
  private def addParameterEnumeration[V](row: Int, pe: Parameter[V]): ComboBox[_] = {
    val label = new Label(pe.name + ":") {
      tooltip = pe.description
    }
    val comboBox = new ComboBox(pe.enabledValues) {
      renderer = ParameterRenderer
    }
    comboBox.selection.item = pe.default
    layout(label) = new Constraints(0, row, 1, 1, 0.0, 0.0, GridBagConstraints.EAST,
      GridBagConstraints.NONE, new Insets(0, 0, 5, 5), 0, 0)
    layout(comboBox) = new Constraints(1, row, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
      GridBagConstraints.HORIZONTAL, new Insets(0, 0, 5, 0), 0, 0)
    comboBox
  }

  def selectedDomain = NumericalDomain.enabledValues(domainComboBox.selection.index)
  
  def getParameters[T <: Target](tgt: T) = {
    val parameters = new Parameters(selectedDomain, tgt)
    parameters.wideningScope = WideningScope(wideningComboBox.selection.index)
    parameters.narrowingStrategy = NarrowingStrategy(narrowingComboBox.selection.index)
    val delay = delayModel.getValue().asInstanceOf[Double].toInt
    if (delay != 0) {
      parameters.wideningFactory = DelayedWideningFactory(DefaultWidening,delay)
    }
    parameters.narrowingFactory = DelayedNarrowingFactory(NoNarrowing,2)
    if (debug.selected) parameters.debugWriter = new java.io.StringWriter
    parameters
  }
}
