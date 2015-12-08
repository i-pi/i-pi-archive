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

__all__=['ForceConstMover']

import numpy as np
import time


from ipi.engine.mover import Mover
from ipi.utils.depend import *
from ipi.utils import units
from ipi.utils.softexit import softexit
from ipi.utils.mintools import min_brent, min_approx, BFGS, L_BFGS, L_BFGS_nls
from ipi.utils.messages import verbosity, warning, info

class ForceConstMover(Mover):
    """Dynamic matrix calculation routine. Computes the force constant matrix and then discrete Fourier 
          interpolation.
    """

    def __init__(self, fixcom=False, fixatoms=None, epsilon=0.001, oldk=0, matrix=np.zeros(0, float)):   
                 
        """Initialises ForceConstMover.
        Args:
        fixcom: An optional boolean which decides whether the centre of mass
                motion will be constrained or not. Defaults to False. 
                THe Hessian is an array of matrix for each beads.        
        """

        super(ForceConstMover,self).__init__(fixcom=fixcom, fixatoms=fixatoms)
      
        #Finite difference option.
        self.epsilon = epsilon
        self.oldk = oldk
	self.matrix = matrix
   
    def bind(self, ens, beads, nm, cell, bforce, bbias, prng):
      
        super(ForceConstMover,self).bind(ens, beads, nm, cell, bforce, bbias, prng)
        if(self.beads.nbeads > 1):
            raise ValueError("Incorrect entry for calculation of force constant matrix")

	if(self.matrix.size  != (beads.q.size * beads.q.size)):
            if(self.matrix.size == 0):
                self.matrix=np.eye(beads.q.size, beads.q.size, 0, float)
            else:
                raise ValueError("Force constant matrix size does not match system size")
            
    def step(self, step=None):
        """Calculates the kth derivative of force by finite differences.            
        """
     
        if(step==None):
            k=0
        else:
            k=step-1

        self.ptime = self.ttime = 0
        self.qtime = -time.time()

        info("\nDynmtarix STEP %d" % step, verbosity.debug)
    
        #initialise des donnes du system compris par IPI
        #if(self.dforces is None):#formations of duplicates
        self.dbeads = self.beads.copy()
        self.dcell = self.cell.copy()
        self.dforces = self.forces.copy(self.dbeads, self.dcell) 
            
        #initialze the vector if doesn't exit or reinitialyze to zero all components a 3N vector
        self.delta = np.zeros(self.beads.nbeads * 3 * self.dbeads.natoms, float)       
        #delta = an array with all ements equal to 0 except that kth element is epsilon.
        self.delta[k] = self.epsilon
        #displaces kth d.o.f by epsilon.                          
        self.dbeads.q = self.beads.q + self.delta  #making it one raw 3N long
        fplus = - depstrip(self.dforces.f).copy()
        #displaces kth d.o.f by -epsilon.      
        self.dbeads.q = self.beads.q - self.delta 
        fminus =  - depstrip(self.dforces.f).copy()
        #computes a row of force-constant matrix
        forces_raw = (fplus-fminus)/(2*self.epsilon * depstrip(self.beads.sm3[-1][k]) * depstrip(self.beads.sm3[-1]) )
        #change the line value or add the line if does not exit to the matrix
        self.matrix[k] = forces_raw
        if (k == 3 * self.dbeads.natoms -1):
            outfile=open('./ForceConstantMatrix.out', 'w+')
            for j in range(0,3 * self.dbeads.natoms):
                print >> outfile, ' '.join(map(str, self.matrix[j]))
            softexit.trigger("Force constant matrix is calculated. Exiting simulation")
