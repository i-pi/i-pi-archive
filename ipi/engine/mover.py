"""Contains the classes that deal with the different dynamics required in
different types of ensembles.

Copyright (C) 2013, Joshua More and Michele Ceriotti

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http.//www.gnu.org/licenses/>.


Holds the algorithms required for normal mode propagators, and the objects to
do the constant temperature and pressure algorithms. Also calculates the
appropriate conserved energy quantity for the ensemble of choice.

"""

__all__=['Mover', 'ReplayMover']

import numpy as np
import time

from ipi.utils.depend import *
from ipi.utils import units
from ipi.utils.softexit import softexit
from ipi.utils.io import read_file
from ipi.utils.io.inputs.io_xml import xml_parse_file
from ipi.utils.units import Constants, unit_to_internal
from ipi.utils.mintools import min_brent, min_approx, BFGS


class Mover(dobject):
   """Base mover calculation class.

   Gives the standard methods and attributes needed in all the
   mover calculation classes.

   Attributes:
      beads: A beads object giving the atoms positions.
      cell: A cell object giving the system box.
      forces: A forces object giving the virial and the forces acting on
         each bead.      
      fixcom: A boolean which decides whether the centre of mass
         motion will be constrained or not.
      fixatoms: A list of atoms that should be held fixed to their 
         initial positions.

   Depend objects:
      none
   """

   def __init__(self, fixcom=False, fixatoms=None):
      """Initialises Mover object.

      Args:
         fixcom: An optional boolean which decides whether the centre of mass
            motion will be constrained or not. Defaults to False.
         fixatoms: A list of atoms that should be held fixed to their 
           initial positions.
      """
            
      self.fixcom = fixcom
      if fixatoms is None: 
         self.fixatoms = np.zeros(0,int)
      else:
         self.fixatoms = fixatoms


   def bind(self, beads, nm, cell, bforce, bbias, prng):
      """Binds beads, cell, bforce, bbias and prng to the calculator.

      This takes a beads object, a cell object, a forcefield object and a
      random number generator object and makes them members of the atom mover caclulator.
      It also then creates the objects that will hold the data needed in the
      ensemble algorithms and the dependency network. Note that the conserved
      quantity is defined in the init, but as each ensemble has a different
      conserved quantity the dependencies are defined in bind.

      Args:
         beads: The beads object from whcih the bead positions are taken.
         nm: A normal modes object used to do the normal modes transformation.
         cell: The cell object from which the system box is taken.
         bforce: The forcefield object from which the force and virial are taken.
         bbias: The forcefield object from which the bias forces are obtained            
         prng: The random number generator object which controls random number
            generation.
      """

      # store local references to the different bits of the simulation
      self.beads = beads
      self.cell = cell
      self.forces = bforce
      self.bias = bbias
      self.prng = prng
      self.nm = nm


   def step(self, step=None):
      """Dummy simulation time step which does nothing."""

      pass      

class ReplayMover(Mover):
   """Calculator object that just loads snapshots from an external file in sequence.

   Has the relevant conserved quantity and normal mode propagator for the
   constant energy ensemble. Note that a temperature of some kind must be
   defined so that the spring potential can be calculated.

   Attributes:
      intraj: The input trajectory file.
      ptime: The time taken in updating the velocities.
      qtime: The time taken in updating the positions.
      ttime: The time taken in applying the thermostat steps.

   Depend objects:
      None really meaningful.
   """

   def __init__(self, fixcom=False, fixatoms=None, intraj=None):
      """Initialises ReplayEnsemble.

      Args:
         dt: The simulation timestep.
         temp: The system temperature.
         fixcom: An optional boolean which decides whether the centre of mass
            motion will be constrained or not. Defaults to False.
         intraj: The input trajectory file.
      """

      super(ReplayMover,self).__init__(fixcom=fixcom, fixatoms=fixatoms)
      if intraj == None:
         raise ValueError("Must provide an initialized InitFile object to read trajectory from")
      self.intraj = intraj
      if intraj.mode == "manual":
         raise ValueError("Replay can only read from PDB or XYZ files -- or a single frame from a CHK file")
      self.rfile = open(self.intraj.value,"r")
      self.rstep = 0

   def step(self, step=None):
      """Does one replay time step."""

      self.ptime = self.ttime = 0
      self.qtime = -time.time()

      
      while True:
       self.rstep += 1
       try:         
         if (self.intraj.mode == "xyz"):            
            for b in self.beads:
               myatoms = read_file("xyz", self.rfile)
               myatoms.q *= unit_to_internal("length",self.intraj.units,1.0)
               b.q[:] = myatoms.q
         elif (self.intraj.mode == "pdb"):
            for b in self.beads:
               myatoms, mycell = read_file("pdb", self.rfile)
               myatoms.q *= unit_to_internal("length",self.intraj.units,1.0)
               mycell.h  *= unit_to_internal("length",self.intraj.units,1.0)
               b.q[:] = myatoms.q
            self.cell.h[:] = mycell.h
         elif (self.intraj.mode == "chk" or self.intraj.mode == "checkpoint"):
            # reads configuration from a checkpoint file
            xmlchk = xml_parse_file(self.rfile) # Parses the file.

            from ipi.inputs.simulation import InputSimulation
            simchk = InputSimulation()
            simchk.parse(xmlchk.fields[0][1])
            mycell = simchk.cell.fetch()
            mybeads = simchk.beads.fetch()
            self.cell.h[:] = mycell.h
            self.beads.q[:] = mybeads.q
            softexit.trigger(" # Read single checkpoint")
       except EOFError:
         softexit.trigger(" # Finished reading re-run trajectory")
       if (step==None or self.rstep>step): break 
      self.qtime += time.time()
