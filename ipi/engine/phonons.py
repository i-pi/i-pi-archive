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

__all__=['DynMatrixMover']

import numpy as np
import time


from ipi.engine.mover import Mover
from ipi.utils.depend import *
from ipi.utils import units
from ipi.utils.softexit import softexit
from ipi.utils.mintools import min_brent, min_approx, BFGS, L_BFGS, L_BFGS_nls
from ipi.utils.messages import verbosity, warning, info

class DynMatrixMover(Mover):
    """Dynamic matrix calculation routine by finite difference.
    """

    def __init__(self, fixcom=False, fixatoms=None, epsilon=0.001, oldk=0, matrix=np.zeros(0, float)):   
                 
        """Initialises DynMatrixMover.
        Args:
        fixcom	: An optional boolean which decides whether the centre of mass
             	  motion will be constrained or not. Defaults to False. 
        matrix	: A 3Nx3N array that stores the dynamic matrix.
        oldk	: An integr that stores the number of rows calculated.
        epsilon: A 3Nx3N array that stores the dynamic matrix.
        """

        super(DynMatrixMover,self).__init__(fixcom=fixcom, fixatoms=fixatoms)
      
        #Finite difference option.
        self.epsilon = epsilon
        self.oldk = oldk
	self.matrix = matrix
        self.sm=np.zeros(0,float)
   
    def bind(self, ens, beads, nm, cell, bforce, bbias, prng):

        #Raises error for nbeads not equal to 1.      
        super(DynMatrixMover,self).bind(ens, beads, nm, cell, bforce, bbias, prng)
        if(self.beads.nbeads > 1):
            raise ValueError("Calculation not possible for number of beads greater than one")

        #Initialises a 3*number of atoms X 3*number of atoms dynamic matrix.
	if(self.matrix.size  != (beads.q.size * beads.q.size)):
            if(self.matrix.size == 0):
                self.matrix=np.eye(beads.q.size, beads.q.size, 0, float)
            else:
                raise ValueError("Force constant matrix size does not match system size")

        #Creates square root mass matrix.
	if(self.sm.size  != (beads.q.size * beads.q.size)):
            if(self.sm.size == 0):
                self.sm=np.diag(beads.m3[-1])
            else:
                raise ValueError("Force constant matrix size does not match system size")

        self.dbeads = self.beads.copy()
        self.dcell = self.cell.copy()
        self.dforces = self.forces.copy(self.dbeads, self.dcell) 
            
    def step(self, step=None):
        """Calculates the kth derivative of force by finite differences.            
        """
     
        if(step==None):
            k=0
        else:
            k=step-1

        self.ptime = self.ttime = 0
        self.qtime = -time.time()

        info("\nDynMatrix STEP %d" % step, verbosity.debug)
    
        #initialise des donnes du system compris par IPI
            
        #initializes the finite deviation
        self.delta = np.zeros(3 * self.beads.natoms, float)       
        self.delta[k] = self.epsilon
        #displaces kth d.o.f by epsilon.                          
        self.dbeads.q = self.beads.q + self.delta  
        plus = - depstrip(self.dforces.f).copy()
        #displaces kth d.o.f by -epsilon.      
        self.dbeads.q = self.beads.q - self.delta 
        minus =  - depstrip(self.dforces.f).copy()
        #computes a row of force-constant matrix
        DynMatrixElement = (plus-minus)/(2*self.epsilon*depstrip(self.beads.sm3[-1][k])*depstrip(self.beads.sm3[-1]))
        #change the line value or add the line if does not exit to the matrix
        self.matrix[k] = DynMatrixElement


        if (k == 3*self.beads.natoms -1):
            eigsys=np.linalg.eig((self.matrix + np.transpose(self.matrix))/2)
            self.q2u=eigsys[0]*self.sm
            self.u2q=np.linalg.inv(self.q2u)
            outfile01=open('./DynMatrix.matrix.out', 'w+')
            outfile02=open('./DynMatrix.eigenvalues.out', 'w+')
            outfile03=open('./DynMatrix.eigenvectors.out', 'w+')
            print >> outfile02, ' '.join(map(str, eigsys[0]))
            for j in range(0,3 * self.dbeads.natoms):
                print >> outfile01, ' '.join(map(str, self.matrix[j]))
                print >> outfile03, ' '.join(map(str, eigsys[1][j]))
            outfile01.close
            outfile02.close
            outfile03.close
            softexit.trigger("Dynamic matrix is calculated. Exiting simulation")
