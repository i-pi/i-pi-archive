import numpy
import math 

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
   
def compute_ih(h) :
   #print "Dummy routine"
   ih = numpy.zeros((3,3), float)
   for i in range(3):
      ih[i,i] = 1.0/h[i,i]
   ih[0,1] = -h[0,1]/(h[0,0] * h[1,1])
   ih[1,2] = -h[1,2]/(h[1,1] * h[2,2])
   ih[0,2] = (h[0,1]*h[1,0] - h[0,2]*h[1,1]) / (h[0,0]*h[1,1]*h[2,2])
   return h
        
class Cell(object):
   """Represents the simulation cell in a periodic system"""

#h is the lattice basis matrix, which will hold the basis vectors in column form
#h is going to be upper triangular, i.e. a_1=(a1,0,0); a_2=(a2x,a2y,0); a_3=(a3x,a3y,a3z)
#p will hold the box momenta, in the same form as for h
#w is probably the barostat mass, or something similar. 
  
   @property
   def h(self): 
      return self.__h
   
   @h.setter
   def h(self, newh): 
      self.__h = newh
      self.__ih = numpy.zeros((3,3),float)
      self.__taint_ih = True

   @property
   def p(self): 
      return self.__p
   
   @p.setter
   def p(self, newp): 
      self.__p=newp

   @property
   def ih(self):       
      if (self.__taint_ih) :
         self.__ih=compute_ih(self.__h)
         self.__taint_ih = False
      return self.__ih
   
   def __init__(self):
      self.__h = numpy.identity(3, float)
      self.__p = numpy.zeros((3,3) ,float)
      self.__taint_ih = True
      self.w = 1.0

   def __str__(self):
      return "    Unit vectors: %s %s %s \n    Momenta: %s %s %s \n    w = %s, volume = %s" % (self.h[0:3,0], self.h[0:3,1], self.h[0:3,2], self.p[0:3,0], self.p[0:3,1], self.p[0:3,2], self.w, self.volume())
      
   def volume(self):
      #return numpy.inner(self.h[0:3,0], numpy.cross(self.h[0:3,1], self.h[0:3,2])) # general matrix
      return self.__h[0,0]*self.__h[1,1]*self.__h[2,2]   # upper-triangular matrix

   def kinetic(self):
      ke = 0.0
      for i in range(3):
         for j in range(3):
            ke += self.__p[i, j]**2
      ke /= 2*self.w
      return ke
      
   def apply_pbc(self, atom):
      s=numpy.dot(self.ih,atom.q)
      for i in range(0,3):
         s[i] = s[i] - round(s[i])
      print s
      atom = numpy.dot(self.h,s)
      print atom
      
   def h2abc(self):
      """
      Returns a description of the cell in terms of the length of the 
      lattice vectors and the angles between them."""
      
      a=self.__h[0,0]; b=math.sqrt(self.__h[0,1]**2+self.__h[1,1]**2);  c=math.sqrt(self.__h[0,2]**2+self.__h[1,2]**2+self.__h[2,2]**2);
      gamma=math.acos(self.__h[0,1]/b); beta=math.acos(self.__h[0,2]/c); alpha=math.acos((self.__h[0,1]*self.__h[0,2])*(self.__h[1,1]*self.__h[1,2])/(b*c));
      return a, b, c, alpha, beta, gamma


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



def print_pdb(atoms, cell):
   a, b, c, alpha, beta, gamma = cell.h2abc()
   alpha *= 180.0/math.pi
   beta  *= 180.0/math.pi
   gamma *= 180.0/math.pi
   
   z = 0 #we need to find out what Z actually is
   print "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%s% 4i" % (a, b, c, alpha, beta, gamma, " P 1       ", z)
   for i in range(0,len(atoms)): 
      print "ATOM  % 5i%4s%c% 3s%c% 4i%c%8.3f%8.3f%8.3f%6.2f%6.2f%2s% 2i" % (i+1, atoms[i].name,'X','  1',' ',1,' ',atoms[i].q[0],atoms[i].q[1],atoms[i].q[2],0.0,0.0,'  ',0)

def read_pdb(filedesc):

#We need a way of using these values to initialise the system

   header = filedesc.readline()
   a = float(header[6:15]);      b = float(header[15:24]);     c = float(header[24:33]);
   alpha = float(header[33:40]); beta = float(header[40:47]);  gamma = float(header[47:54]);

   body = filedesc.readline()
   while body != '':
      name = body[12:16]
      x = float(body[30:38])
      y = float(body[38:46])
      z = float(body[46:54])
      pos = numpy.array([x, y, z])
      body = filedesc.readline()
   atoms = ''
   cell = ''
   return atoms, cell
