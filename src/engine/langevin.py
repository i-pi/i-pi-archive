import numpy, math, random
import thermostat
from utils.depend import *

class Thermo_Langevin(thermostat.Thermostat):     
   
   def compute_T(self):
      return math.exp(-self.dt/self.tau)
      
   def compute_S(self):      
      return math.sqrt(self.k_Boltz*self.temp*(1-self.T**2))
   
   def __init__(self, temp = 1.0, dt = 1.0, tau = 1.0, econs=0.0):
      super(Thermo_Langevin,self).__init__(temp,dt,econs)
      self.tau=depend(value=tau,name='tau')
      self.T=depend(name='T',func=self.compute_T)
      self.S=depend(name='S',func=self.compute_S)      
      
      self.getdesc('tau').add_dependant(self.getdesc('T'))
      self.getdesc('dt').add_dependant(self.getdesc('T'))
      self.getdesc('T').add_dependant(self.getdesc('S'))
      self.getdesc('temp').add_dependant(self.getdesc('S'))
                  
      self.k_Boltz = 1.0

   def step(self, atom):
      sm=math.sqrt(atom.mass)
      self.econs+=atom.kin
      for i in range(3):
         atom.p[i] = self.T*atom.p[i] + sm*self.S*random.gauss(0.0, 1.0)
      self.econs-=atom.kin      

   def cell_step(self, cell):
      for i in range(3):
         for j in range(i,3):
            cell.p[i,j] = self.T*cell.p[i,j] + self.S*math.sqrt(cell.w)*random.gauss(0.0, 1.0)

