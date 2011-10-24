import numpy, math, random
from thermostat import *
from utils.depend import *
from utils import units

class langevin(thermostat):     
   def compute_T(self):
      return math.exp(-self.dt.get()/self.tau.get())
      
   def compute_S(self):      
      return math.sqrt(units.kb*self.temp.get()*(1-self.T.get()**2))
   
   def compute_smass(self):
      sm=numpy.zeros(3*len(self.atoms))
      for i in range(len(self.atoms)):
         sm[3*i]=sm[3*i+1]=sm[3*i+2]=math.sqrt(self.atoms[i].mass.get())
      return sm

   def compute_sw(self):
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
      # stores links to the momentum array, and constructs a sqrt(mass) dependency array
      #TODO makes sure that the shapes of p and atoms are compatible
      self.p=p
      self.atoms=atoms
      self.cell=cell
      
      self.smass=depend(name='smass',func=self.compute_smass)
      for at in self.atoms: at.mass.add_dependant(self.smass)
      
      self.sw=depend(name='sw',func=self.compute_sw, deplist=[self.cell.w])
      

   def step(self):
      sm=self.smass.get();  p=self.p.get();  T=self.T.get(); S=self.S.get();
      econs=self.econs.get()
      for i in range(len(p)):
         p[i]/=sm[i]
         econs+=p[i]*p[i]*0.5
         p[i] = T*p[i]+ S*random.gauss(0.0,1.0);
         econs-=p[i]*p[i]*0.5
         p[i]*=sm[i]
      self.econs.set(econs)
      self.p.taint(taintme=False)

   def cell_step(self):
      #TODO define properly within the new framework
      sw = self.sw.get(); T=self.T.get(); S=self.S.get();   p=self.cell.p.get()
      self.econs.set(self.econs.get()+self.cell.kin.get())
      for i in range(3):
         for j in range(i,3): 
            p[i,j] = T*p[i,j] + S*sw*random.gauss(0.0, 1.0)
      self.cell.pc.set(self.cell.pc.get()*T+S*sw*random.gauss(0.0, 1.0))
      self.cell.p.taint(taintme=False)           
      self.econs.set(self.econs.get()-self.cell.kin.get())      

