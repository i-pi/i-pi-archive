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
"""

__all__=['Dynmatrix']

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
      oldj   : value of j in previous step.
      oldk   : value of k in previous step.
   """

   def __init__(self, fixcom=False, fixatoms=None, epsilon=0.001,oldj=0,oldk=0,noldbead=0,oldhessian=np.zeros(0, float)) :   
                 
      """Initialises Dynmatrix.
      Args:
         fixcom: An optional boolean which decides whether the centre of mass
            motion will be constrained or not. Defaults to False. 
            THe Hessian is an array of matrix for each beads.        
      """

      super(Dynmatrix,self).__init__(fixcom=fixcom, fixatoms=fixatoms)
      
      #Finite difference option.
      self.epsilon = epsilon
      self.oldj = oldj
      self.oldk = oldk
      self.hessian = oldhessian
      self.noldbead =noldbead
   
	def bind(self, ens, beads, nm, cell, bforce, bbias, prng):
      
		super(Dynmatrix,self).bind(ens, beads, nm, cell, bforce, bbias, prng)
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
            
	def step(self, k, step=None):
      """Calculates the kth derivative of force by finite differences.            
      """
      
		self.ptime = self.ttime = 0
		self.qtime = -time.time()

		info("\nDynmtarix STEP %d" % step, verbosity.debug)
	
	#initialise des donnes du system compris par IPI
		if(self.dforces is None) :#formations of duplicates
			self.dbeads = self.beads.copy()
			self.dcell = self.cell.copy()
			self.dforces = self.bforce.copy(self.dbeads, self.dcell) 
	#delta = an array with all ements equal to 0 except that kth element is epsilon.
    #displaces kth d.o.f by epsilon.	  	  			  
		self.dbeads.q = self.beads.q + gen_delta(k)  
		fplus = - destrip(self.dforces.f)[k]
	# displaces kth d.o.f by -epsilon.	  
		self.dbeads.q = self.beads.q - gen_delta(k) 
		fminus =  - destrip(self.dforces.f)[k]
	# computes the derivative			  
		forces_raw = (fplus-fminus)/(2*self.epsilon)
	#change the array value or add the array if does not exit
		if self.hessian.size < k:
			self.hessian = np.column_stack([sef.hessian, forces_raw])
        else:
			self.hessian[:,k] = forces_raw
		  
	#update the row k per bead
	   
	def gen_delta(self, k):
		"""Change the delta which becomes an array with 1 to the kth index
		and zero elsewhere
		"""
	#initialze the vector if doesn't exit or reinitialyze to zero all components a 3N vector
		#self.delta = np.zeros(beads.q.size, float)
		self.delta = np.zeros(self.dbeads.nbeads * 3 * self.dbeads.natoms, float)	   
		self.delta[k]=self.epsilon
		
	
		
	

	
	   
		 
	   
			  




