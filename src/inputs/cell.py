"""Deals with creating the cell class.

Generates an cell class from a cell vector.

Classes:
   InputCell: Deals with creating the Cell object from a file, and
      writing the checkpoints.
"""

import numpy as np
import math
from copy import copy
import utils.mathtools as mt
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

   attribs = copy(InputArray.attribs)
   attribs["mode"] = (InputAttribute, { "dtype"  : str,
                                        "default": "h",
                                        "options": ["h", "abc", "abcABC"],
                                        "help"   : "This decides whether the system box is created from a cell parameter matrix, or from the side lengths and angles between them. If 'mode' is 'h', then 'cell' takes a 3*3 cell vector matrix. If 'mode' is 'abcABC', then 'cell' takes an array of 6 floats, the first three being the length of the sides of the system parallelopiped, and the last three being the angles (in degrees) between those sides. Angle A corresponds to the angle between sides b and c, and so on for B and C. If mode is 'abc', then this is the same as ffor 'abcABC', but the cell is assumed to be orthorhombic."} )

   default_help = "Deals with the cell parameters. Takes as arguments either the cell vector matrix, or the unit cell side lengths and the angles between them."
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
      self.mode.store("h")  # we always store the cell matrix

   def fetch(self):
      """Creates a cell object.

      Returns:
         A cell object of the appropriate type and with the appropriate
         properties given the attributes of the InputCell object.
      """

      h = super(InputCell,self).fetch()
      h.shape = (3,3)

      return Cell(h=h)

   def check(self):
      """Checks that the array which has been read represent a valid cell
      matrix.
      """

      super(InputCell,self).check()

      h = self.value
      if self.mode.fetch() == "h":
         if h.size != 9:
            raise ValueError("Cell objects must contain a 3x3 matrix describing the cell vectors.")
      elif self.mode.fetch() == "abc":
         if h.size != 3:
            raise ValueError("If you are initializing cell from cell side lengths you must pass the 'cell' tag an array of 3 floats.")
         else:
            h = mt.abc2h(h[0], h[1], h[2], math.pi/2, math.pi/2, math.pi/2)
            super(InputCell,self).store(h)
            self.shape.store((3,3))
      elif self.mode.fetch() == "abcABC":
         if h.size != 6:
            raise ValueError("If you are initializing cell from cell side lengths and angles you must pass the 'cell' tag an array of 6 floats.")
         else:
            h = mt.abc2h(h[0], h[1], h[2], h[3]*math.pi/180.0, h[4]*math.pi/180.0, h[5]*math.pi/180.0)
            super(InputCell,self).store(h)
            self.shape.store((3,3))

      h.shape = (9,)
      if not (h[3] == 0.0 and h[6] == 0.0 and h[7] == 0.0):
         print "Warning: cell vector matrix must be upper triangular, all elements below the diagonal being set to zero."
         h[3] = 0
         h[6] = 0
         h[7] = 0
