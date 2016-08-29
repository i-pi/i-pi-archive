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


from ipi.engine.motion import Motion
from ipi.utils.depend import *
from ipi.utils import units
from ipi.utils.softexit import softexit
from ipi.utils.messages import verbosity, warning, info

class DynMatrixMover(Motion):
    """Dynamic matrix calculation routine by finite difference.
    """

    def __init__(self, fixcom=False, fixatoms=None, mode='std', energy_shift=0.0, pos_shift=0.001, output_shift=0.000,
                 dynmat=np.zeros(0, float), dynmat_r=np.zeros(0, float), prefix="", asr="none"):   
                 
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
        if self.mode == "fd":
            self.phononator = self.FDPhononator()
        elif self.mode == "posref":
            self.phononator = self.FDPhononator()
        elif self.mode == "enrgref":
            self.phononator = self.FDPhononator()

        self.deltaw = output_shift
        self.deltax = pos_shift
        self.deltae = energy_shift 
        self.dynmatrix = dynmat
        self.dynmatrix_r = dynmat_r
        self.frefine = False
        self.U = None
        self.V = None
        self.prefix = prefix
        if self.prefix == "":
            self.prefix = "PHONONS"
        self.asr = asr
   
    def bind(self, ens, beads, nm, cell, bforce, prng):

        super(DynMatrixMover,self).bind(ens, beads, nm, cell, bforce, prng)

        #Raises error for nbeads not equal to 1.    
        if(self.beads.nbeads > 1):
            raise ValueError("Calculation not possible for number of beads greater than one")

        #Initialises a 3*number of atoms X 3*number of atoms dynamic matrix.
        if(self.dynmatrix.size  != (beads.q.size * beads.q.size)):
            if(self.dynmatrix.size == 0):
                if(self.mode == "fd"):
                    self.dynmatrix=np.zeros((beads.q.size, beads.q.size), float)
                    self.dynmatrix_r=np.zeros((beads.q.size, beads.q.size), float)
                else:
                    raise ValueError("Force constant matrix size not found")
            else:
                raise ValueError("Force constant matrix size does not match system size")
 
        self.phononator.bind(self)
        self.ism = 1/np.sqrt(depstrip(self.beads.m3[-1]))
        self.m = depstrip(self.beads.m)
    
    def step(self, step=None):
        """Computes one step of phonon computation. """
        if (step < 3*self.beads.natoms):
            self.phononator.step(step)
        else:
            print "ASR", self.asr
	    rdyn = self.apply_asr()
            self.printall(self.prefix, rdyn)
            softexit.trigger("Dynamic matrix is calculated. Exiting simulation")                    

    class DummyPhononator(dobject):
        """ No-op phononator """

        def __init__(self):
            pass
 
        def bind(self, dm):
            """ Reference all the variables for simpler access."""
 
            self.beads = dm.beads
            self.cell = dm.cell
            self.ensemble = dm.ensemble
            self.forces = dm.forces
            self.deltaw = dm.deltaw
            self.deltax = dm.deltax
            self.deltae = dm.deltae
            self.dynmatrix = dm.dynmatrix

            self.dbeads = self.beads.copy()
            self.dcell = self.cell.copy()
            self.dforces = self.forces.copy(self.dbeads, self.dcell)
            self.ism = 1/np.sqrt(depstrip(self.beads.m3[-1]))
            self.m = depstrip(self.beads.m)
 
        def step(self, step=None):
            """Dummy simulation time step which does nothing."""
            pass

    class FDPhononator(DummyPhononator):
        """ Finite dinnerence phonon evaluator.
        """

        def step(self, step=None):
            """Computes one row of the dynamic matrix by finite difference."""
 
            #initializes the finite deviation
            dev = np.zeros(3 * self.beads.natoms, float)
            dev[step] = self.deltax
            #displaces kth d.o.f by delta.                          
            self.dbeads.q = self.beads.q + dev
            plus = - depstrip(self.dforces.f).copy()
            #displaces kth d.o.f by -delta.      
            self.dbeads.q = self.beads.q - dev
            minus =  - depstrip(self.dforces.f).copy()
            #computes a row of force-constant matrix
            dmrow = (plus-minus)/(2*self.deltax)*self.ism[step]*self.ism
            self.dynmatrix[step] = dmrow

    def printall(self, prefix, dmatx, deltaw=0.0):
        """ Prints out diagnostics for a given dynamical matrix. """

        dmatx = dmatx + np.eye(len(dmatx))*deltaw
        if deltaw != 0.0 :
            wstr = " !! Shifted by %e !!" % (deltaw)
        else:
            wstr = ""
        # prints out the dynamical matrix
        outfile=open(prefix+'.dynmat', 'w')
        print >> outfile, "# Dynamical matrix (atomic units)"+wstr
        for i in range(3 * self.beads.natoms):
            print >> outfile, ' '.join(map(str, dmatx[i]))
        outfile.close()
        
        # also prints out the Hessian
        outfile=open(prefix+'.hess', 'w')
        print >> outfile, "# Hessian matrix (atomic units)"+wstr
        for i in range(3 * self.beads.natoms):
            print >> outfile, ' '.join(map(str, dmatx[i]/(self.ism[i]*self.ism)))
        outfile.close()
        
        eigsys=np.linalg.eigh(dmatx)        
        # prints eigenvalues & eigenvectors
        outfile=open(prefix+'.eigval', 'w') 
        print >> outfile, "# Eigenvalues (atomic units)"+wstr
        print >> outfile, '\n'.join(map(str, eigsys[0]))
        outfile.close()
        outfile=open(prefix+'.eigvec', 'w')        
        print >> outfile, "# Eigenvector  matrix (normalized)"
        for i in range(0,3 * self.beads.natoms):
            print >> outfile, ' '.join(map(str, eigsys[1][i]))
        outfile.close()
        
        eigmode = 1.0*eigsys[1]
        for i in range(0,3 * self.beads.natoms):
            eigmode[i] *= self.ism[i]
        for i in range(0,3 * self.beads.natoms):
            eigmode[:,i]/=np.sqrt(np.dot(eigmode[:,i],eigmode[:,i]))
        outfile=open(prefix+'.mode', 'w')        
        print >> outfile, "# Phonon modes (mass-scaled)"
        for i in range(0,3 * self.beads.natoms):
            print >> outfile, ' '.join(map(str, eigmode[i]))
        outfile.close()
            
    def sstep(self, step=None):
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
                if edelta > 100*self.deltax:  edelta= 100*self.deltax          
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
            #computes a row of the refined dynmatrix, in the basis of the eigenvectors of the first dynmatrix            
            
            dmrowk = (plus-minus)/(2*edelta/vknorm)
            
            self.dynmatrix_r[k] = np.dot(self.V.T, dmrowk)
                        
        if k >= 3*self.beads.natoms-1:
            # symmetrize and apply chosen Acoustic Sum Rule
            rdyn = self.asr_apply(self.dynmatrix)
            
            if not self.frefine:            
                self.printall(self.prefix, rdyn, self.deltaw)
                if self.mode=="std":
                    softexit.trigger("Dynamic matrix is calculated. Exiting simulation")                    
            else:
                rdyn = self.asr_apply(self.dynmatrix_r)
                
                self.printall(self.prefix + "-R", rdyn, self.deltaw)
                # transform in Cartesian basis
                rdyn = np.dot(self.U,np.dot(rdyn,np.transpose(self.U)))
                self.printall(self.prefix + "-RC", rdyn, self.deltaw)
                softexit.trigger("Dynamic matrix is calculated. Exiting simulation")                    

    def apply_asr(self):
        """
        Removes the translations and/or rotations depending on the asr mode.
        """
        if(self.asr=="none"):
            pass

        if(self.asr=="crystal"):
            #Computes the centre of mass.
            com=np.dot(np.transpose(self.beads.q.reshape((self.beads.natoms,3))),self.m)/self.m.sum()
            qminuscom= self.beads.q.reshape((self.beads.natoms,3)) - com
            #Computes the moment of inertia tensor.
            moi=np.zeros((3,3), float)
            for k in range(self.beads.natoms):
                moi =moi -np.dot(np.cross(qminuscom[k],np.identity(3)),np.cross(qminuscom[k],np.identity(3)))*self.m[k]

            U=(np.linalg.eig(moi))[1]
            R=np.dot(qminuscom,U)
            D=np.zeros((3,3*self.beads.natoms),float)

            #Computes the vectors along rotations.
            D[0]=np.tile([1,0,0],self.beads.natoms)/self.ism
            D[1]=np.tile([0,1,0],self.beads.natoms)/self.ism
            D[2]=np.tile([0,0,1],self.beads.natoms)/self.ism

            #Computes unit vecs.
            for k in range(3):
                D[k]=D[k]/np.linalg.norm(D[k])

            #Computes the transformation matrix.
            transfmatrix = np.eye(3*self.beads.natoms)-np.dot(D.T,D)
            return np.dot(transfmatrix.T,np.dot(self.dynmatrix,transfmatrix))

        elif(self.asr=="poly"):
            #Computes the centre of mass.
            com=np.dot(np.transpose(self.beads.q.reshape((self.beads.natoms,3))),self.m)/self.m.sum()
            qminuscom= self.beads.q.reshape((self.beads.natoms,3)) - com
            #Computes the moment of inertia tensor.
            moi=np.zeros((3,3), float)
            for k in range(self.beads.natoms):
                moi =moi -np.dot(np.cross(qminuscom[k],np.identity(3)),np.cross(qminuscom[k],np.identity(3)))*self.m[k]

            U=(np.linalg.eig(moi))[1]
            R=np.dot(qminuscom,U)
            D=np.zeros((6,3*self.beads.natoms),float)

            #Computes the vectors along translations and rotations.
            D[0]=np.tile([1,0,0],self.beads.natoms)/self.ism
            D[1]=np.tile([0,1,0],self.beads.natoms)/self.ism
            D[2]=np.tile([0,0,1],self.beads.natoms)/self.ism
            for i in range(3*self.beads.natoms):
                iatom=i/3
                idof=np.mod(i,3)
                D[3,i]=(R[iatom,1]*U[idof,2] - R[iatom,2]*U[idof,1])/self.ism[i]
                D[4,i]=(R[iatom,2]*U[idof,0] - R[iatom,0]*U[idof,2])/self.ism[i]
                D[5,i]=(R[iatom,0]*U[idof,1] - R[iatom,1]*U[idof,0])/self.ism[i]

            #Computes unit vecs.
            for k in range(6):
                D[k]=D[k]/np.linalg.norm(D[k])

            #Computes the transformation matrix.
            transfmatrix = np.eye(3*self.beads.natoms)-np.dot(D.T,D)
            return np.dot(transfmatrix.T,np.dot(self.dynmatrix,transfmatrix))
