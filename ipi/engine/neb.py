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

__all__=['NEBMover']

import numpy as np
import time


from ipi.engine.mover import Mover
from ipi.utils.depend import *
from ipi.utils import units
from ipi.utils.softexit import softexit
from ipi.utils.io import read_file
from ipi.utils.io.inputs.io_xml import xml_parse_file
from ipi.utils.units import Constants, unit_to_internal
from ipi.utils.mintools import min_brent, min_approx, BFGS

class NEBBFGSMover(object):
   """ Creation of the multi-dimensional function that will be minimized"""

   def __init__(self):
      self.x0 = self.d = self.xold = None

   def bind(self, ens):
      self.dbeads = ens.beads.copy()
      self.dcell = ens.cell.copy()
      self.dforces = ens.forces.copy(self.dbeads, self.dcell)

   def __call__(self, x):
       
      self.dbeads.q = x
      bq = depstrip(self.dbeads.q).copy()
      bf = depstrip(self.dforces.f).copy()
      
      nimg = self.dbeads.nbeads
      nat = self.dbeads.natoms
      
      # get tangents
      btau = np.zeros((nimg, 3*nat), float)
      for ii in range(1,nimg-1):
          d1 = bq[ii]-bq[ii-1]
          d2 = bq[ii+1]-bq[ii]
          btau[ii]= d1/np.linalg.norm(d1)+d2/np.linalg.norm(d2)
          btau[ii] *= 1.0/np.linalg.norm(btau)
        
      # get perpendicular forces 
      for ii in range(1,nimg-1):
          bf[ii] = bf[ii] - np.dot(bf[ii],btau[ii]) * btau[ii]
          print np.linalg.norm(bf[ii])
          
      # adds the spring forces           
      for ii in range(1,nimg-1):
          print np.dot( btau[ii], ( bq[ii+1]+bq[ii-1]-2*bq[ii] )  ) 
          bf[ii] += self.neb_kappa * btau[ii]*np.dot( btau[ii], ( bq[ii+1]+bq[ii-1]-2*bq[ii] )  ) 

 
      e = 0.0
      g = -bf
      return e, g  
      
class NEBMover(Mover):
   """Nudged elastic band routine.

   Attributes:

   """

   def __init__(self, fixcom=False, fixatoms=None,
             mode="sd", 
             grad_tolerance=1.0e-6, maximum_step=100.0,
             cg_old_force=np.zeros(0, float),
             cg_old_direction=np.zeros(0, float),
             invhessian=np.eye(0), 
             ls_options={ "tolerance": 1e-5,  "iter": 100.0 , "step": 1e-3, "adaptive":1.0 } ,
             tolerances = {"energy" : 1e-5, "force": 1e-5, "position": 1e-5}
             ) :   
                 
      """Initialises GeopMover.

      Args:
         fixcom: An optional boolean which decides whether the centre of mass
            motion will be constrained or not. Defaults to False.         
      """

      super(NEBMover,self).__init__(fixcom=fixcom, fixatoms=fixatoms)
      
      # optimization options
      self.ls_options = ls_options
      self.tolerances = tolerances
      self.mode=mode
      self.grad_tol = grad_tolerance
      self.max_step = maximum_step
      self.cg_old_f = cg_old_force
      self.cg_old_d = cg_old_direction
      self.invhessian = invhessian
      self.neb_kappa = 1e0 # TODO make this a parameter!
      self.bfgs = NEBBFGSMover()
         
   
   def bind(self, beads, nm, cell, bforce, bbias, prng):
      
      super(NEBMover,self).bind(beads, nm, cell, bforce, bbias, prng)
      if self.cg_old_f.size != beads.q.size :
         if self.cg_old_f.size == 0: 
            self.cg_old_f = np.zeros(beads.q.size, float)
         else: 
            raise ValueError("Conjugate gradient force size does not match system size")
      if self.cg_old_d.size != beads.q.size :
         if self.cg_old_d.size == 0: 
            self.cg_old_d = np.zeros(beads.q.size, float)
         else: 
            raise ValueError("Conjugate gradient direction size does not match system size")
      
      self.bfgs.bind(self)
      
   def step(self, step=None):
      """Does one simulation time step."""

      self.ptime = self.ttime = 0
      self.qtime = -time.time()

      bq = depstrip(self.beads.q).copy()
      bf = depstrip(self.forces.f).copy()
      
      nimg = self.beads.nbeads
      nat = self.beads.natoms
      
      # get tangents
      btau = np.zeros((nimg, 3*nat), float)
      for ii in range(1,nimg-1):
          d1 = bq[ii]-bq[ii-1]
          d2 = bq[ii+1]-bq[ii]
          btau[ii]= d1/np.linalg.norm(d1)+d2/np.linalg.norm(d2)
          btau[ii] *= 1.0/np.linalg.norm(btau)
        
      # get perpendicular forces 
      for ii in range(1,nimg-1):
          bf[ii] = bf[ii] - np.dot(bf[ii],btau[ii]) * btau[ii]
          print np.linalg.norm(bf[ii])
          
      # adds the spring forces           
      for ii in range(1,nimg-1):
          print np.dot( btau[ii], ( bq[ii+1]+bq[ii-1]-2*bq[ii] )  ) 
          bf[ii] += self.neb_kappa * btau[ii]*np.dot( btau[ii], ( bq[ii+1]+bq[ii-1]-2*bq[ii] )  ) 

      self.beads.q += bf *0.05
      
      self.qtime += time.time()
