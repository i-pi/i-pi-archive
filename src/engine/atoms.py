import numpy

class Atom:
    """Represent an atom, with position, velocity, mass and related properties"""
    
#qpslice holds position and momentum data. 
#qslice holds a reference to the position data, pslice to the momentum data.

    def __init__(self, qpslice, name="X", mass=1.0):
        self.q = qpslice[0:3,0]
        self.p = qpslice[0:3,1]
        self.mass = mass
	self.name = name

    def __str__(self):
        return "    Name = %s, q = %s, p = %s " % (self.name, self.q, self.p)

    def pot(self):
	pot = 0.0
        return pot

    def kinetic(self):
	ke = 0.0
	for i in range(3):
	    ke += self.p[i]**2 / (2.0*self.mass)
	return ke
        
class Cell:
   """Represents the simulation cell in a periodic system"""

#h is the lattice basis matrix, which will hold the basis vectors in column form
#p will hold the box momenta, in the same form as for h
#w is probably the barostat mass, or something similar. 

   def __init__(self):
      self.h = numpy.identity(3, float)
      self.p = numpy.zeros((3,3) ,float)
      self.w = 1.0

   def __str__(self):
      return "    Unit vectors: %s %s %s \n    Momenta: %s %s %s \n    w = %s, volume = %s" % (self.h[0:3,0], self.h[0:3,1], self.h[0:3,2], self.p[0:3,0], self.p[0:3,1], self.p[0:3,2], self.w, self.volume())
      
   def volume(self):
      return numpy.inner(self.h[0:3,0], numpy.cross(self.h[0:3,1], self.h[0:3,2]))

   def kinetic(self):
      ke = 0.0
      for i in range(3):
         for j in range(3):
            ke += self.p[i, j]**2
      ke /= 2*self.w
      return ke
        
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
