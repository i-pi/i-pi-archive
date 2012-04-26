"""Deals with creating the cell class.

Generates an cell class either from a set of positions and momenta, or from 
a configuration file.

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
      
class InputCell(Input):
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
      file_units: An optional string giving the length units that the file is
         specified by. Defaults to ''.
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
            "file_units": (InputValue, {"dtype"    : str,
                                        "default"  : "",
                                        "help"     : "The units in which the lengths in the configuration file are given.",
                                        "options"  : [unit for unit in UnitMap["length"]] }),
            "from_file": (InputValue, {"dtype"     : str,
                                       "default"   : "",
                                       "help"      : "A file from which to take the cell parameters from."}) }
   attribs={ "flexible" : (InputValue, {"dtype"    : bool, 
                                        "default"  : False,
                                        "help"     : "Whether the cell parameters can change during the simulation."}) }

   default_help = "Deals with the cell parameters, and stores their momenta in flexible cell calculations."
   default_label = "CELL"
    
   def store(self, cell, filename=""):
      """Takes a Cell instance and stores of minimal representation of it.

      Args:
         cell: A cell object.
         filename: An optional float giving a file to read the cell dimensions
            from. Defaults to ''.
      """

      super(InputCell,self).store(cell)
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
         properties given the attributes of the InputCell object.
      """

      super(InputCell,self).fetch()
      if (self.flexible.fetch()): 
         cell = CellFlexi(h=self.h.fetch(), m=self.m.fetch())
         if not self.h0._explicit:
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

      super(InputCell,self).check()
      if self.from_file._explicit:
         
         filename=self.from_file.fetch()
         ext=filename[len(filename)-3:]
         if (ext == "pdb"):
            myatoms, mycell = utils.io.io_pdb.read_pdb(open(self.from_file.fetch(),"r"))
            mycell.h *= UnitMap["length"][self.file_units.fetch()]
            mycell.h0 *= UnitMap["length"][self.file_units.fetch()]
         else:
            raise ValueError("Unrecognized extension for cell configuration")

         if not self.h._explicit:
            self.h.store(mycell.h)
         if not self.h0._explicit:
            self.h0.store(mycell.h0)
         if not self.m._explicit:
            self.m.store(mycell.m)
         self.from_file.store("")

      h = self.h.fetch()
      h0 = self.h0.fetch()
      if not (h.shape == (3,3) and h0.shape == (3,3)):
         raise ValueError("Incorrect shape for cell vector matrices.")
      if not (h[1,0] == h[2,0] == h[2,1] == h0[1,0] == h0[2,0] == h0[2,1] == 0):
         print "Warning: cell vector matrix must be upper triangular, all elements below the diagonal being set to zero."
         h[1,0] = 0
         h[2,0] = 0
         h[2,1] = 0
         h0[1,0] = 0
         h0[2,0] = 0
         h0[2,1] = 0
