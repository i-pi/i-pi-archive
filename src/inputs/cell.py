"""Deals with creating the cell class.

Generates an cell class from a cell vector.

Classes:
   InputCell: Deals with creating the Cell object from a file, and
      writing the checkpoints.
"""

import numpy as np
import utils.io.io_pdb, utils.io.io_xyz
from engine.cell import *
from utils.inputvalue import *
from utils.units import UnitMap

__all__ = [ 'InputCell' ]

class InputCell(InputArray):
   """Cell input class.

   Handles generating the appropriate cell class from the xml input file,
   and generating the xml checkpoint tags and data from an instance of the
   object.
   """

   default_help = "Deals with the cell parameters."
   default_label = "CELL"

   def __init__(self, help=None, dimension=None, units=None, default=None, dtype=None):
      """ Initializes an InputCell object by just calling the parent
          with appropriate arguments. """

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

      h=super(InputCell,self).fetch();
      h.shape=(3,3)

      return Cell(h=h)

   def check(self):
      """Checks that the array which has been read represent a valid cell
      matrix.
      """

      super(InputCell,self).check()

      h = self.value
      if h.size != 9:
         raise ValueError("Cell objects must contain a 3x3 matrix describing the cell vectors.")

      h.shape = (9,)
      if not (h[3] == 0.0 and h[6] == 0.0 and h[7] == 0.0):
         print "Warning: cell vector matrix must be upper triangular, all elements below the diagonal being set to zero."
         h[3] = 0
         h[6] = 0
         h[7] = 0
