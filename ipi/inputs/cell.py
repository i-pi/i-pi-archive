"""Deals with creating the cell class.

Copyright (C) 2013, Joshua More and Michele Ceriotti

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


Generates an cell class from a cell vector.

Classes:
   InputCell: Deals with creating the Cell object from a file, and
      writing the checkpoints.
"""

import numpy as np
from copy import copy
from ipi.engine.cell import *
from ipi.utils.inputvalue import *
from ipi.utils.units import UnitMap
from ipi.utils.messages import verbosity, warning

__all__ = [ 'InputCell' ]

class InputCell(InputArray):
   """Cell input class.

   Handles generating the appropriate cell class from the xml input file,
   and generating the xml checkpoint tags and data from an instance of the
   object.
   """

   attribs = copy(InputArray.attribs)

   default_help = "Deals with the cell parameters. Takes as array which can be used to initialize the cell vector matrix."
   default_label = "CELL"

   def __init__(self, help=None, dimension=None, units=None, default=None, dtype=None):
      """Initializes InputCell.

      Just calls the parent initialization function with appropriate arguments.
      """

      super(InputCell,self).__init__(dtype=float, dimension="length", default=default, help=help)

   def store(self, cell):
      """Takes a Cell instance and stores of minimal representation of it.

      Args:
         cell: A cell object.
      """

      super(InputCell,self).store(cell.h)
      self.shape.store((3,3))

   def fetch(self):
      """Creates a cell object.

      Returns:
         A cell object of the appropriate type and with the appropriate
         properties given the attributes of the InputCell object.
      """

      h = super(InputCell,self).fetch()
      h.shape = (3,3)

      return Cell(h=h)
