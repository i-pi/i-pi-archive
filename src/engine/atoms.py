import numpy

class Atom:
    """Represent an atom, with position, velocity, mass and related properties"""
    
    def __init__(self, qpslice, name="X", mass=1.0):
        self.q = qpslice[0:3,0]
        self.p = qpslice[0:3,1]
        self.mass = mass;  self.name = name
    def __str__(self):
        return "%4s %s %s " % (self.name, self.q, self.p)
        
class Cell:
   """Represents the simulation cell in a periodic system"""
   def __init__(self):
      h = numpy.zeros((3,3) ,float)
      p = numpy.zeros((3,3) ,float)
      w = 1.0
      
        
class System:
    """
    Represents a simulation cell. 
    Includes the cell parameters, the atoms and the like. """
    
    def __init__(self, natoms=1):
        self.natoms=natoms
        self.__qp=numpy.zeros((3*natoms,2),float) 
        self.__qp[:,0]=numpy.arange(0,3*natoms)
        self.q=self.__qp[:,0]
        self.__qp[:,1]=numpy.arange(0,3*natoms)*0.01
        self.p=self.__qp[:,1]
        self.atoms = [ Atom(self.__qp[3*i:3*(i+1),0:2]) for i in range(natoms) ]

    def __str__(self):
        rstr="ATOMS ("+str(self.natoms)+")\n"
        for i in range(0,self.natoms): 
            rstr=rstr+str(self.atoms[i])+"\n"
        return rstr
        
    def step(self,dt):
        self.q+=self.p*dt
   
