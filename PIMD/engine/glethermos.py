import numpy, math, random
from thermostat import *
from utils.depend import *
from utils import units
from utils.mathtools import *

class glethermo(thermostat):

   def compute_T(self):
      return matrix_exp(-self.dt.get()*self.A.get())
      
   def compute_S(self):
      S=units.kb*self.C
      S-=numpy.dot(self.T.get(),numpy.dot(S,numpy.transpose(self.T.get())))
      return S
   
   def compute_smass(self):
      sm=numpy.zeros(3*len(self.atoms))
      for i in range(len(self.atoms)):
         sm[3*i]=sm[3*i+1]=sm[3*i+2]=math.sqrt(self.atoms[i].mass.get())
      return sm

   def compute_sw(self):
      return math.sqrt(self.cell.w.get())
     
   def get_ns(self): return self.A.get().shape[1]; 

   def __init__(self, temp = 1.0, dt = 1.0, A=numpy.zeros(1,float), C=Null, econs=0.0):
      super(glethermo,self).__init__(temp,dt,econs)
      
      self.A=depend(value=numpy.copy(A),name='A')
      self.ns=depend(name='ns',func=self.get_ns)

      # temp may be used to initialize C, but is then ignored
      if (C is Null) : LC=temp*numpy.identity(self.ns.get(),float)
      else: LC=numpy.copy(C)
      self.C=depend(value=LC,name='C')

      self.T=depend(name='T',func=self.compute_T)
      self.S=depend(name='S',func=self.compute_S)      
      
      self.A.add_dependant(self.ns)
      self.A.add_dependant(self.T)
      self.dt.add_dependant(self.T)
      self.T.add_dependant(self.S)
      self.C.add_dependant(self.S)

   def bind(self, atoms, p, cell):
      # stores links to the momentum array, and constructs a sqrt(mass) dependency array
      #TODO makes sure that the shapes of p and atoms are compatible
      self.p=p
      self.atoms=atoms
      self.cell=cell
      
      self.smass=depend(name='smass',func=self.compute_smass)
      for at in self.atoms: at.mass.add_dependant(self.smass)
      
      self.sw=depend(name='sw',func=self.compute_sw, deplist=[self.cell.w])
      self._ps=numpy.zeros((self.ns.get()+1,len(self.atoms)*3), float)
      self._xi=numpy.zeros((self.ns.get()+1,len(self.atoms)*3), float)
      #TODO here I should also initialize the s's?

   def step(self):
      sm=self.smass.get();  p=self.p.get();  T=self.T.get(); S=self.S.get();

      econs=self.econs.get()

      p[:]/=sm[:]

      econs+=numpy.inner(p,p)*0.5
      self._ps[1,:]=p[:]   # gets mass-scaled p inside the first row of ps   
      
      for el in self._xi.flat: el=random.gauss(0.0,1.0)

      self._ps[:]=numpy.dot(T,self._ps)+numpy.dot(S,self.xi)

      p[:]=self._ps[1,:]
      econs-=numpy.inner(p,p)*0.5
      
      p[:]*=sm[:]

      self.econs.set(econs)
      self.p.taint(taintme=False)

