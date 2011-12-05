import numpy, math, random
from thermostat import *
from utils.depend import *
from utils import units

class PILE(thermostat):     
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

      gamma = self.gamma.get_array()
      T = numpy.zeros(len(gamma))
      for i in range(len(gamma)):
         T[i] = math.exp(-self.dt.get()*gamma[i])
      
      return T
      
   def compute_S(self):      
      """Calculates S in p(0) = T*p(dt) + S*random.gauss()"""

      T = self.T.get_array()
      S = numpy.zeros(len(T))
      for i in range(len(T)):
         S[i] = math.sqrt(units.kb*self.temp.get()*(1-T[i]**2))

      return S
   
   def compute_smass(self):
      """Calculates sqrt(mass)"""

      atoms = self.syst.atoms
      sm=numpy.zeros(len(atoms))
      for i in range(len(atoms)):
         sm[i]=math.sqrt(atoms[i].mass.get())
      return sm

   def compute_sw(self):
      """Calculates sqrt(cell mass)"""

      return math.sqrt(self.syst.cell.w.get())

   def compute_gamma(self):
      
      gamma_vec = 2*self.syst.n_frequencies.get_array()
      gamma_vec[0] = 1.0/self.tau_0.get()
      
      return gamma_vec
     
   def __init__(self, temp = 1.0, dt = 1.0, tau_0 = 1.0, econs=0.0):
      super(PILE,self).__init__(temp,dt,econs)

      self.tau_0 = depend(value=tau_0, name='tau_0')
      
   def bind(self, system):
      """Binds the appropriate system objects to the thermostat, such that
         the thermostat step automatically updates the same velocities of the 
         atoms and the cell"""
      # stores links to the momentum array, and constructs a sqrt(mass) dependency array
      #TODO makes sure that the shapes of p and atoms are compatible

      self.syst = system
      
      self.smass=depend(value=numpy.zeros(len(system.atoms)), name='smass',func=self.compute_smass)
      for at in self.syst.atoms:
         at.mass.add_dependant(self.smass)
      
      self.sw=depend(name='sw',func=self.compute_sw, deplist=[self.syst.cell.w])
      
      self.gamma=depend(value=numpy.zeros(len(system.systems)), func=self.compute_gamma, name='gamma', deplist=[self.tau_0, self.syst.n_frequencies])
      self.T=depend(value=numpy.zeros(len(system.systems)), name='T',func=self.compute_T)
      self.S=depend(value=numpy.zeros(len(system.systems)), name='S',func=self.compute_S)      
      
      self.gamma.add_dependant(self.T)
      self.dt.add_dependant(self.T)
      self.T.add_dependant(self.S)
      self.temp.add_dependant(self.S)      

   def step(self):
      """Updates the atom velocities with a langevin thermostat"""

      sm=self.smass.get_array(); T=self.T.get_array(); S=self.S.get_array();
      econs=self.econs.get()
      p_tilde = numpy.dot(self.syst.trans_mat.get_array(), self.syst.p.get_array())

      for i in range(len(sm)):
         p_tilde_i = numpy.array(p_tilde[:,3*i:3*(i+1)])
         p_tilde_i/=sm[i]
         for j in range(len(self.syst.systems)):
            for k in range(3):
               econs += 0.5*p_tilde_i[j,k]**2
         for j in range(len(self.syst.systems)):
            p_tilde_i[j,:] = T[j]*p_tilde_i[j,:] + S[j]*random.gauss(0.0,1.0);
         for j in range(len(self.syst.systems)):
            for k in range(3):
               econs -= 0.5*p_tilde_i[j,k]**2
         p_tilde_i*=sm[i]
         p_tilde[:,3*i:3*(i+1)] = p_tilde_i

      self.syst.p.get_array()[:] = numpy.dot(numpy.transpose(self.syst.trans_mat.get_array()), p_tilde)
      self.econs.set(econs)
      self.syst.p.taint(taintme=False)

#   def NST_cell_step(self):
#      """Updates the cell velocities with a langevin thermostat in the
#         NST ensemble"""
#
#      sw = self.sw.get(); T=self.T.get(); S=self.S.get(); p=numpy.array(self.cell.p.get(), ndmin=3)
#      self.econs.set(self.econs.get()+self.cell.kin.get())
#      for i in range(3):
#         for j in range(i,3): 
#            for k in range(len(p[:,i,j])):
#               p[k,i,j] = T*p[k,i,j] + S*sw*random.gauss(0.0, 1.0)
#      p.shape=self.cell.p.get().shape
#      self.cell.p.get()[:] = p
#      self.cell.p.taint(taintme=False)           
#      self.econs.set(self.econs.get()-self.cell.kin.get())      
#
#   def NPT_cell_step(self):
#      """Updates the cell velocities with a langevin thermostat in the 
#         NPT ensemble"""
##TODO make this work for a Cell_vec as well as a Cell
#
#      sw = self.sw.get(); T=self.T.get(); S=self.S.get()
#      self.econs.set(self.econs.get()+self.cell.kin.get())
#      self.cell.pc.set(self.cell.pc.get()*T+S*sw*random.gauss(0.0, 1.0))
#      self.econs.set(self.econs.get()-self.cell.kin.get())      
