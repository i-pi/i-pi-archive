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

Algorithms implemented by Michele Ceriotti and Benjamin Helfrecht, 2015

"""

__all__=['GeopMover']

import numpy as np
import time


from ipi.engine.mover import Mover
from ipi.utils.depend import *
from ipi.utils import units
from ipi.utils.softexit import softexit
from ipi.utils.mintools import min_brent, min_approx, BFGS, L_BFGS, L_BFGS_nls
from ipi.utils.messages import verbosity, warning, info

class Dynmatrix(Mover):
   """Dynamic matrix calculation routine. Computes the force constant matrix and then discrete Fourier 
      interpolation.

   Attributes:
      epsilon: list of previous gradients (g_n+1 - g_n) for L-BFGS. Number of entries = corrections      
      oldj   : value of j in previous step.
      oldk   : value of k in previous step.
   """

   def __init__(self, fixcom=False, fixatoms=None, epsilon=0.01,oldj=0,oldk=0,oldhessian=np.zeros(0,float) :   
                 
      """Initialises GeopMover.

      Args:
         fixcom: An optional boolean which decides whether the centre of mass
            motion will be constrained or not. Defaults to False.         
      """

      super(GeopMover,self).__init__(fixcom=fixcom, fixatoms=fixatoms)
      
      #Finite difference option.
      self.epsilon = epsilon
      self.oldj = oldj
      self.oldk = oldk
      self.hessian = oldhessian
   
   def bind(self, ens, beads, nm, cell, bforce, bbias, prng):
      
      super(GeopMover,self).bind(ens, beads, nm, cell, bforce, bbias, prng)
      #if self.cg_old_f.shape != beads.q.size :
      #   if self.cg_old_f.size == 0: 
      #      self.cg_old_f = np.zeros(beads.q.size, float)
      #   else: 
      #      raise ValueError("Conjugate gradient force size does not match system size")
      #if self.cg_old_d.size != beads.q.size :
      #   if self.cg_old_d.size == 0: 
      #      self.cg_old_d = np.zeros(beads.q.size, float)
      #   else: 
      #      raise ValueError("Conjugate gradient direction size does not match system size")
            
   def step(self, j, k, step=None):
      """Calculates the jth derivative of force on the kth atom."""

      
      self.ptime = self.ttime = 0
      self.qtime = -time.time()

      info("\nMD STEP %d" % step, verbosity.debug)

      if(self.dforces is None) :
         self.dbeads = self.beads.copy()
         self.dcell = self.cell.copy()
         self.dforces = self.bforce.copy(self.dbeads, self.dcell) 

      epsilon = 1.0e-3

      #compute element of Hessian.
      #self.dbeads.q = self.beads.q + epsilon# move forward (should hardcode or input displacement)
      #fplus = depstrip(self.dforces.f).copy()
      #self.dbeads.q = self.beads.q - epsilon # move forward (should hardcode or input displacement)
      #fminus = depstrip(self.dforces.f).copy()
      #hessian = 2*(fminus - fplus)/2.0/epsilon

      #append it to self.Hessian



