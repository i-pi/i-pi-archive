import numpy
import math
from io_system import *
from atoms import *
from cell import *

class System:
    """
    Represents a simulation cell. 
    Includes the cell parameters, the atoms and the like. """

#__qp holds all the positions and momenta for all the atoms in the simulation
#q and p hold the positions and momenta, respectively.
#P_ext will be the external load.
#we will probably have to redo the initialisation step, so that it makes sense physically.
#step will eventually call the forces from the external program and then do the propagation step. At the moment we simply take free particle trajectories, to test the theory.
    
    def __init__(self, natoms = 1, temp = 5.0):
        self.natoms=natoms
	self.temp = temp
        self.__qp=numpy.zeros((3*natoms,2),float) 
        self.__qp[:,0]=numpy.arange(0,3*natoms)
        self.q=self.__qp[:,0]
        self.__qp[:,1]=numpy.arange(0,3*natoms)*0.01
        self.p=self.__qp[:,1]
        self.atoms = [ Atom(self.__qp[3*i:3*(i+1),0:2]) for i in range(natoms) ] #Creates a list of atoms from the __qp array
	self.cell = Cell()
	self.P_ext = numpy.zeros(3,float)

    def __str__(self):
        rstr="ATOMS ("+str(self.natoms)+"):\n"
        for i in range(0,self.natoms): 
            rstr=rstr+str(self.atoms[i])+"\n"
	rstr = rstr + "Cell:\n" + str(self.cell)
        return rstr
        
    def step(self,dt):
        self.q+=self.p*dt

    def kinetic(self):
	ke = 0.0
	for i in range(self.natoms):
	    ke += self.atoms[i].kinetic()
	ke += self.cell.kinetic()
	return ke
