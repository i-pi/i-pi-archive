import numpy, math, random
from thermostat import *
from utils.depend import *
from utils import units

class langevin(thermostat):     
   """Represent a langevin thermostat for constant T simulations.
      Contains: temp = temperature, dt = time step, econs =
      change in the kinetic energy due to the thermostat,
      tau = thermostat mass, sm = sqrt(mass), (T,S) = thermostat parameters
      Initialised by: thermo = langevin(temp, dt, tau, econs)
      temp = temperature, default = 1.0
      dt = time step, default = 1.0
      tau = thermostat mass, default = 1.0
      econs = conserved energy quantity, default = 0.0"""

   def compute_T(self):
      """Calculates T in p(0) = T*p(dt) + S*random.gauss()"""

      return math.exp(-self.dt.get()/self.tau.get())
      
   def compute_S(self):      
      """Calculates S in p(0) = T*p(dt) + S*random.gauss()"""

      return math.sqrt(units.kb*self.temp.get()*(1-self.T.get()**2))
   
   def compute_smass(self):
      """Calculates sqrt(mass)"""

      sm=numpy.zeros(3*len(self.atoms))
      for i in range(len(self.atoms)):
         sm[3*i]=sm[3*i+1]=sm[3*i+2]=math.sqrt(self.atoms[i].mass.get())
      return sm

   def compute_sw(self):
      """Calculates sqrt(cell mass)"""

      return math.sqrt(self.cell.w.get())
     
  
   def __init__(self, temp = 1.0, dt = 1.0, tau = 1.0, econs=0.0):
      super(langevin,self).__init__(temp,dt,econs)
      
      self.tau=depend(value=tau,name='tau')
      self.T=depend(name='T',func=self.compute_T)
      self.S=depend(name='S',func=self.compute_S)      
      
      self.tau.add_dependant(self.T)
      self.dt.add_dependant(self.T)
      self.T.add_dependant(self.S)
      self.temp.add_dependant(self.S)      

   def bind(self, atoms, p, cell):
      """Binds the appropriate system objects to the thermostat, such that
         the thermostat step automatically updates the same velocities of the 
         atoms and the cell"""
      # stores links to the momentum array, and constructs a sqrt(mass) dependency array
      #TODO makes sure that the shapes of p and atoms are compatible

      self.p=p
      self.atoms=atoms
      self.cell=cell
      
      self.smass=depend(name='smass',func=self.compute_smass)
      for at in self.atoms: at.mass.add_dependant(self.smass)
      
      self.sw=depend(name='sw',func=self.compute_sw, deplist=[self.cell.w])
      

   def step(self):
      """Updates the atom velocities with a langevin thermostat"""

      sm=self.smass.get(); p=numpy.array(self.p.get(), ndmin=2); T=self.T.get(); S=self.S.get();
      econs=self.econs.get()
      for i in range(len(sm)):
         p[:,i]/=sm[i]
         econs+=numpy.dot(p[:,i],p[:,i])*0.5
         for j in range(len(p[:,i])):
            p[j,i] = T*p[j,i]+ S*random.gauss(0.0,1.0);
         econs-=numpy.dot(p[:,i],p[:,i])*0.5
         p[:,i]*=sm[i]

      p.shape=self.p.get().shape
      self.p.get()[:] = p
      self.econs.set(econs)
      self.p.taint(taintme=False)

   def NST_cell_step(self):
      """Updates the cell velocities with a langevin thermostat in the
         NST ensemble"""

      sw = self.sw.get(); T=self.T.get(); S=self.S.get(); p=numpy.array(self.cell.p.get(), ndmin=3)
      self.econs.set(self.econs.get()+self.cell.kin.get())
      for i in range(3):
         for j in range(i,3): 
            for k in range(len(p[:,i,j])):
               p[k,i,j] = T*p[k,i,j] + S*sw*random.gauss(0.0, 1.0)
      p.shape=self.cell.p.get().shape
      self.cell.p.get()[:] = p
      self.cell.p.taint(taintme=False)           
      self.econs.set(self.econs.get()-self.cell.kin.get())      

   def NPT_cell_step(self):
      """Updates the cell velocities with a langevin thermostat in the 
         NPT ensemble"""
#TODO make this work for a Cell_vec as well as a Cell

      sw = self.sw.get(); T=self.T.get(); S=self.S.get()
      self.econs.set(self.econs.get()+self.cell.kin.get())
      self.cell.pc.set(self.cell.pc.get()*T+S*sw*random.gauss(0.0, 1.0))
      self.econs.set(self.econs.get()-self.cell.kin.get())      
