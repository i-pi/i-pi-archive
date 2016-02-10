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


from ipi.engine.motion.motion import Motion
from ipi.utils.depend import *
from ipi.utils import units
from ipi.utils.softexit import softexit
from ipi.utils.mintools import min_brent, min_approx, BFGS, L_BFGS, L_BFGS_nls
from ipi.utils.messages import verbosity, warning, info

class DynMatrixMover(Motion):
    """Dynamic matrix calculation routine by finite difference.
    """

    def __init__(self, fixcom=False, fixatoms=None, mode='std', energy_shift=1.0, pos_shift=0.001, 
                 dynmat=np.zeros(0, float), dynmat_r=np.zeros(0, float)):   
                 
        """Initialises DynMatrixMover.
        Args:
        fixcom	: An optional boolean which decides whether the centre of mass
             	  motion will be constrained or not. Defaults to False. 
        matrix	: A 3Nx3N array that stores the dynamic matrix.
        oldk	: An integr that stores the number of rows calculated.
        delta: A 3Nx3N array that stores the dynamic matrix.
        """

        super(DynMatrixMover,self).__init__(fixcom=fixcom, fixatoms=fixatoms)
      
        #Finite difference option.
        self.mode = mode
        self.deltax = pos_shift
        self.deltae = energy_shift        
        self.dynmatrix = dynmat
        self.dynmatrix_r = dynmat_r
        self.frefine = False
        self.U = None
        self.V = None
   
    def bind(self, ens, beads, nm, cell, bforce, prng):

        super(DynMatrixMover,self).bind(ens, beads, nm, cell, bforce, prng)

        #Raises error for nbeads not equal to 1.    
        if(self.beads.nbeads > 1):
            raise ValueError("Calculation not possible for number of beads greater than one")

        #Initialises a 3*number of atoms X 3*number of atoms dynamic matrix.
        if(self.dynmatrix.size  != (beads.q.size * beads.q.size)):
            if(self.dynmatrix.size == 0):
                self.dynmatrix=np.zeros((beads.q.size, beads.q.size), float)
                self.dynmatrix_r=np.zeros((beads.q.size, beads.q.size), float)
            else:
                raise ValueError("Force constant matrix size does not match system size")

        self.dbeads = self.beads.copy()
        self.dcell = self.cell.copy()
        self.dforces = self.forces.copy(self.dbeads, self.dcell)         
        self.ism = 1/np.sqrt(depstrip(beads.m3[-1]))        
    
    def printall(self, prefix, dmatx):
        """ Prints out diagnostics for a given dynamical matrix. """

        # prints out the dynamical matrix
        outfile=open(prefix+'.dynmat', 'w')
        print >> outfile, "# Dynamical matrix (atomic units)"
        for i in range(3 * self.dbeads.natoms):
            print >> outfile, ' '.join(map(str, dmatx[i]))
        outfile.close()
        
        # also prints out the Hessian
        outfile=open(prefix+'.hess', 'w')
        print >> outfile, "# Hessian matrix (atomic units)"
        for i in range(3 * self.dbeads.natoms):
            print >> outfile, ' '.join(map(str, dmatx[i]/(self.ism[i]*self.ism)))
        outfile.close()
        
        eigsys=np.linalg.eigh(dmatx)        
        # prints eigenvalues & eigenvectors
        outfile=open(prefix+'.eigval', 'w') 
        print >> outfile, "# Eigenvalues (atomic units)"
        print >> outfile, '\n'.join(map(str, eigsys[0]))
        outfile.close()
        outfile=open(prefix+'.eigvec', 'w')        
        print >> outfile, "# Eigenvector  matrix (normalized)"
        for i in range(0,3 * self.dbeads.natoms):
            print >> outfile, ' '.join(map(str, eigsys[1][i]))
        outfile.close()
        
        eigmode = 1.0*eigsys[1]
        for i in range(0,3 * self.dbeads.natoms):
            eigmode[i] *= self.ism[i]
        for i in range(0,3 * self.dbeads.natoms):
            eigmode[:,i]/=np.sqrt(np.dot(eigmode[:,i],eigmode[:,i]))
        outfile=open(prefix+'.mode', 'w')        
        print >> outfile, "# Phonon modes (mass-scaled)"
        for i in range(0,3 * self.dbeads.natoms):
            print >> outfile, ' '.join(map(str, eigmode[i]))
        outfile.close()
        
        
                    
            
        
    def step(self, step=None):
        """Calculates the kth derivative of force by finite differences.            
        """
     
        if(step == None):
            k = 0
        elif(step < 3*self.beads.natoms): # round one
            k = step            
        else: # round two
            k = step-3*self.beads.natoms
            self.frefine = True
        print "K ", k
        self.ptime = self.ttime = 0
        self.qtime = -time.time()
        info("\nDynMatrix STEP %d" % step, verbosity.debug)
        
        dev = np.zeros(3 * self.beads.natoms, float)       
        if not self.frefine: # fill up the density matrix
            #initializes the finite deviation
            dev[:] = 0
            dev[k] = self.deltax

            #displaces kth d.o.f by delta.                          
            self.dbeads.q = self.beads.q + dev  
            plus = - depstrip(self.dforces.f).copy()

            #displaces kth d.o.f by -delta.      
            self.dbeads.q = self.beads.q - dev 
            minus =  - depstrip(self.dforces.f).copy()

            #computes a row of force-constant matrix
            dmrow = (plus-minus)/(2*self.deltax)*self.ism[k]*self.ism
            self.dynmatrix[k] = dmrow
        else:
            # if needed, computes the eigenvalues of the base matrix
            if self.U is None:
                eigsys=np.linalg.eigh(self.dynmatrix)        
                self.w2 = eigsys[0]
                self.U = eigsys[1]
                self.V = eigsys[1].copy()
                for i in xrange(len(self.V)): self.V[:,i]*=self.ism
                print "U", self.U
                print "V", self.V

            #initializes the finite deviation along one of the (mass-scaled) eigenvectors
            vknorm = np.sqrt(np.dot(self.V[:,k],self.V[:,k]))
            dev = np.real(self.V[:,k]/vknorm)
            
            if self.mode=="nrg":
                edelta = vknorm*np.sqrt(self.deltae*2.0/abs(self.w2[k]))                
            else:
                edelta = self.deltax
            dev *= edelta
            print "displace by", edelta

            #displaces by -delta along kth normal mode.
            self.dbeads.q = self.beads.q + dev
            plus = - depstrip(self.dforces.f).copy().flatten()
            #displaces by -delta along kth normal mode.
            self.dbeads.q = self.beads.q - dev
            minus =  - depstrip(self.dforces.f).copy().flatten()
            #computes a row of the refin    ed dynmatrix, in the basis of the eigenvectors of the first dynmatrix            
            
            dmrowk = (plus-minus)/(2*edelta/vknorm)
            
            self.dynmatrix_r[k] = np.dot(self.V.T, dmrowk)
                        
        if k >= 3*self.beads.natoms-1:
            # symmetrize
            self.dynmatrix = 0.5*(self.dynmatrix+np.transpose(self.dynmatrix))
            
            if not self.frefine:            
                self.printall("PHONONS", self.dynmatrix)
                if self.mode=="std":
                    softexit.trigger("Dynamic matrix is calculated. Exiting simulation")                    
            else:
                self.dynmatrix_r = 0.5*(self.dynmatrix_r+np.transpose(self.dynmatrix_r))
                self.printall("PHONONS-R", self.dynmatrix_r)
                # transform in Cartesian basis
                self.dynmatrix_r = np.dot(self.U,np.dot(self.dynmatrix_r,np.transpose(self.U)))
                self.printall("PHONONS-RC", self.dynmatrix_r)
                softexit.trigger("Dynamic matrix is calculated. Exiting simulation")                    
                
