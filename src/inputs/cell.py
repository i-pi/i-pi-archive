"""Deals with creating the cell class.

Generates an cell class either from a set of positions and momenta, or from 
a configuration file.

Classes:
   RestartCell: Deals with creating the Cell object from a file, and 
      writing the checkpoints.
"""

import numpy as np
from utils.io.io_pdb import *
from engine.cell import *
from utils.inputvalue import *

__all__ = [ 'RestartCell' ]
      
class RestartCell(Input):
   """Cell input class.

   Handles generating the appropriate cell class from the xml input file,
   and generating the xml checkpoint tags and data from an instance of the 
   object.

   Attributes:
      m: An optional float giving the mass of the cell. Defaults to 0.0.
      h: An optional array giving the system box. Defaults to a 3*3 identity 
         matrix.
      h0: An optional array giving the reference box. Defaults to an array 
         of zeros.
      p: An optional array giving the lattice momentum matrix. Defaults to an 
         array of zeros.
      P: An optional float giving the conjugate momentum to the volume 
         fluctuations. Defaults to 0.0.
      init_temp: An optional float to give the effective temperature that the 
         cell velocities should be initialised to. Defaults to -1.0.
      from_file: An optional string giving the name of a pdb file containing
         the initial cell and atom positions. Defaults to ''.
      flexible: A boolean giving whether the cell will be allowed to change 
         shape. Defaults to False.
   """

   fields={ "m" : (InputValue, {"dtype"      : float, 
                                  "default"    : 0.0,
                                  "help"       : "The 'mass' of the cell, used in constant pressure simulations.",
                                  "dimension"  : "mass"}),
            "h" : (InputArray, {"dtype"      : float,
                                  "default"    : np.zeros((3,3)),
                                  "help"       : "The cell vector matrix",
                                  "dimension"  : "length"}), 
            "h0" : (InputArray, {"dtype"     : float,
                                   "default"   : np.zeros((3,3)), 
                                   "help"      : "The reference cell vector matrix. Defined as the unstressed equilibrium cell.",
                                   "dimension" : "length"}),
            "p" : (InputArray, {"dtype"      : float,
                                  "default"    : np.zeros((3,3),float),
                                  "help"       : "The cell 'momenta' matrix, used in constant pressure simulations.",
                                  "dimension"  : "momentum"}),
            "P" : (InputValue, {"dtype"      : float,
                                  "default"    : 0.0,
                                  "help"       : "The scalar cell 'momentum', used in constant pressure simulations.",
                                  "dimension"  : "momentum"}),
            "init_temp": (InputValue, {"dtype"     : float, 
                                         "default"   : -1.0,
                                         "help"      : "The temperature at which the initial velocity distribution is taken, if applicable.",
                                         "dimension" : "temperature"}),
            "from_file": (InputValue, {"dtype"     : str,
                                         "default"   : "",
                                         "help"      : "A file from which to take the cell parameters from.",}) }
   attribs={ "flexible" : (InputValue, {"dtype"    : bool, 
                                          "default"  : False,
                                          "help"     : "Whether the cell parameters can change during the simulation."}) }
    
   def __init__(self, cell=None):
      """Initialises RestartCell.

      Args:
         cell: A Cell object from which to initialise from.
      """

      super(RestartCell,self).__init__()
      if not cell is None:
         self.store(cell)
      
   def store(self, cell, filename=""):
      """Takes a Cell instance and stores of minimal representation of it.

      Args:
         cell: A cell object.
         filename: An optional float giving a file to read the cell dimensions
            from. Defaults to ''.
      """

      super(RestartCell,self).store(cell)
      self.from_file.store(filename)
      self.flexible.store(type(cell) is CellFlexi)
      self.m.store(cell.m)
      self.h.store(cell.h)
      self.h0.store(cell.h0)
      if type(cell) is CellFlexi:
         self.p.store(cell.p)
      else:
         self.P.store(cell.P[0])
      self.from_file.store(cell.from_file)
      
   def fetch(self):
      """Creates a cell object.

      Returns: 
         A cell object of the appropriate type and with the appropriate 
         properties given the attributes of the RestartCell object.
      """

      super(RestartCell,self).fetch()
      if (self.flexible.fetch()): 
         cell = CellFlexi(h=self.h.fetch(), m=self.m.fetch())
         if det_ut3x3(self.h0.fetch()) == 0.0:
            cell.h0 = cell.h
         cell.p = self.p.fetch()
      else:
         cell = CellRigid(h=self.h.fetch(), m=self.m.fetch())
         cell.P = self.P.fetch()
      cell.from_file = self.from_file.fetch()
      return cell
      
   def check(self):
      """Function that deals with optional arguments.
      
      Deals with the init_temp and from_file arguments, and uses
      them to initialise some of the cell parameters depending on which ones
      have been specified explicitly.
      """

      super(RestartCell,self).check()
      if self.from_file.fetch() != "":
         myatoms, mycell = utils.io.io_pdb.read_pdb(open(self.from_file.fetch(),"r")) 
         if (self.h0.fetch() == np.zeros((3,3))).all():
            self.h.store(mycell.h)
         if (self.h0.fetch() == np.zeros((3,3))).all():
            self.h0.store(mycell.h0)
         if self.m.fetch() == 0.0:
            self.m.store(mycell.m)
         self.from_file.store("")
