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

__all__=['Mover', 'ReplayMover', 'GeopMover', 'GeoMin']

import numpy as np
import time

from ipi.utils.depend import *
from ipi.utils import units
from ipi.utils.softexit import softexit
from ipi.utils.io.io_xyz import read_xyz
from ipi.utils.io.io_pdb import read_pdb
from ipi.utils.io.io_xml import xml_parse_file
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
               myatoms = read_xyz(self.rfile)
               myatoms.q *= unit_to_internal("length",self.intraj.units,1.0)
               b.q[:] = myatoms.q
         elif (self.intraj.mode == "pdb"):
            for b in self.beads:
               myatoms, mycell = read_pdb(self.rfile)
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


class LineMover(object):
   """Creation of the one-dimensional function that will be minimized"""
   
   def __init__(self):
      self.x0 = self.d = None

   def bind(self, ens):
      self.dbeads = ens.beads.copy()
      self.dcell = ens.cell.copy()
      self.dforces = ens.forces.copy(self.dbeads, self.dcell)      
      
   def set_dir(self, x0, mdir):      
      self.x0 = x0.copy()
      self.d = mdir.copy()/np.sqrt(np.dot(mdir.flatten(),mdir.flatten()))
      if self.x0.shape != self.d.shape: raise ValueError("Incompatible shape of initial value and displacement direction")      
   
   def __call__(self, x):
            
      self.dbeads.q = self.x0 + self.d * x
      e = self.dforces.pot
      g = - np.dot(depstrip(self.dforces.f).flatten(),self.d.flatten())
      return e, g
   
class BFGSMover(object):
   """ Creation of the multi-dimensional function that will be minimized"""

   def __init__(self):
      self.x0 = self.d = None

   def bind(self, ens):
      self.dbeads = ens.beads.copy()
      self.dcell = ens.cell.copy()
      self.dforces = ens.forces.copy(self.dbeads, self.dcell)

   def __call__(self, x):
      self.dbeads.q = x
      e = self.dforces.pot
      g = - self.dforces.f
      return e, g

class GeoMin(object):
    """ Stores options and utilities for geometry optimization. """
    
    def __init__(self, mode="sd", lin_iter=1000, lin_step=1.0e-3, lin_tol=1.0e-6, grad_tol=1.0e-6, lin_auto=True, max_step=100.0,
             cg_old_f=np.zeros(0, float),
             cg_old_d=np.zeros(0, float),
             invhessian=np.eye(0)):
        self.mode=mode
        self.lin_tol = lin_tol
        self.grad_tol = grad_tol
        self.lin_iter = lin_iter
        self.max_step = max_step
        self.lin_step = lin_step
        self.lin_auto = lin_auto
        self.cg_old_f = cg_old_f
        self.cg_old_d = cg_old_d
        self.invhessian = invhessian

class GeopMover(Mover):
   """Geometry optimization routine. Will start with a dumb steepest descent,
   and then see to include more and more features - getting at least to BFGS.

   Attributes:

   """

   def __init__(self, fixcom=False, fixatoms=None, geop = None):      
      """Initialises GeopMover.

      Args:
         fixcom: An optional boolean which decides whether the centre of mass
            motion will be constrained or not. Defaults to False.         
      """

      super(GeopMover,self).__init__(fixcom=fixcom, fixatoms=fixatoms)
      self.lm = LineMover()
      self.bfgsm = BFGSMover()
      if geop is None: geop = GeoMin()
      self.mo = geop
   
   def bind(self, beads, nm, cell, bforce, bbias, prng):
      
      super(GeopMover,self).bind(beads, nm, cell, bforce, bbias, prng)
      if self.mo.cg_old_f.size != beads.q.size :
         if self.mo.cg_old_f.size == 0: 
            self.mo.cg_old_f = np.zeros(beads.q.size, float)
         else: 
            raise ValueError("Conjugate gradient force size does not match system size")
      if self.mo.cg_old_d.size != beads.q.size :
         if self.mo.cg_old_d.size == 0: 
            self.mo.cg_old_d = np.zeros(beads.q.size, float)
         else: 
            raise ValueError("Conjugate gradient direction size does not match system size")
            
      self.lm.bind(self)
      self.bfgsm.bind(self)
      
   def step(self, step=None):
      """Does one simulation time step."""

      self.ptime = self.ttime = 0
      self.qtime = -time.time()

      print "\nMD STEP %d\n" % step

      if (self.mo.mode == "bfgs"):

          # BFGS Minimization
          # Initialize approximate Hessian inverse and direction
          # to the steepest descent direction
          if step == 0:# or np.sqrt(np.dot(self.bfgsm.d, self.bfgsm.d)) == 0.0:
              self.bfgsm.d = depstrip(self.forces.f) / np.sqrt(np.dot(self.forces.f.flatten(), self.forces.f.flatten()))
              self.mo.invhessian = np.eye(len(self.beads.q.flatten()))
              print "invhessian:", self.mo.invhessian
              #invhessian = np.eye(n)

          # Current energy, forces, and function definitions
          # for use in BFGS algorithm

          u0, du0 = (self.forces.pot, - self.forces.f)
          print "u0:", u0
          print "du0:", du0
          
          # Do one iteration of BFGS, return new point, function value,
          # derivative value, and current Hessian to build upon
          # next iteration
          print "BEFORE"
          print "self.beads.q:", self.beads.q
          print "self.bfgsm.d:", self.bfgsm.d
          print "invhessian:", self.mo.invhessian
          self.beads.q, fx, self.bfgsm.d, self.mo.invhessian = BFGS(self.beads.q, self.bfgsm.d, self.bfgsm, fdf0=(u0, du0), invhessian=self.mo.invhessian, max_step=self.mo.max_step, tol=self.mo.lin_tol, grad_tol=self.mo.grad_tol, itmax=self.mo.lin_iter)  #TODO: make object for inverse hessian and direction if necessary
          print "AFTER"
          print "self.beads.q", self.beads.q
          print "self.bfgsm.d", self.bfgsm.d
          print "invhessian", self.mo.invhessian

      # Routine for steepest descent and conjugate gradient
      else:
          if (self.mo.mode == "sd" or step == 0): 
          
              # Steepest descent minimization
              # gradf1 = force at current atom position
              # dq1 = direction of steepest descent
              # dq1_unit = unit vector of dq1
              gradf1 = dq1 = depstrip(self.forces.f)

              # move direction for steepest descent and 1st conjugate gradient step
              dq1_unit = dq1 / np.sqrt(np.dot(gradf1.flatten(), gradf1.flatten())) 
      
          else:
          
              # Conjugate gradient, Polak-Ribiere
              # gradf1: force at current atom position
              # gradf0: force at previous atom position
              # dq1 = direction to move
              # dq0 = previous direction
              # dq1_unit = unit vector of dq1
              gradf0 = self.mo.cg_old_f
              dq0 = self.mo.cg_old_d
              gradf1 = depstrip(self.forces.f)
              beta = np.dot((gradf1.flatten() - gradf0.flatten()), gradf1.flatten()) / (np.dot(gradf0.flatten(), gradf0.flatten()))
              dq1 = gradf1 + max(0.0, beta) * dq0
              dq1_unit = dq1 / np.sqrt(np.dot(dq1.flatten(), dq1.flatten()))

          self.mo.cg_old_d[:] = dq1    # store force and direction for next CG step
          self.mo.cg_old_f[:] = gradf1
   
          if (len(self.fixatoms)>0):
              for dqb in dq1_unit:
                  dqb[self.fixatoms*3] = 0.0
                  dqb[self.fixatoms*3+1] = 0.0
                  dqb[self.fixatoms*3+2] = 0.0
      
          self.lm.set_dir(depstrip(self.beads.q), dq1_unit)

          # reuse initial value since we have energy and forces already
          u0, du0 = (self.forces.pot, np.dot(depstrip(self.forces.f.flatten()), dq1_unit.flatten()))

          (x, fx) = min_brent(self.lm, fdf0=(u0, du0), x0=0.0, tol=self.mo.lin_tol, itmax=self.mo.lin_iter, init_step=self.mo.lin_step) 

          if self.mo.lin_auto: self.mo.lin_step = x # automatically adapt the search step for the next iteration
      
          self.beads.q += dq1_unit * x


      self.qtime += time.time()
