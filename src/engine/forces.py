import math, numpy
import io_system
from utils.depend import *

class forcefield(object):
   """Creates an interface between the force, potential and virial calculation
      and the wrapper.
      Contains: _pot = potential, _f = forces, _vir = virial, _ufv = 
      (potential, forces, virial)
      Initialised by: ffield = forcefield()"""

   _pot=0.0
   _f=numpy.zeros(0,float)
   _vir=numpy.zeros((3,3),float)
   
   def __init__(self, pars=dict()):
      self._ufv=(self._pot, self._f, self._vir)

   def bind(self, cell, atoms, pot, f, vir):
      """Binds the appropriate system objects to the force calculation, such
         that the potential, force and virial calculation automatically updates
         the same quantities in the system"""

      self.cell=cell
      self.atoms=atoms      # expects a shallow copy of the system's atoms list            
      self._f.resize(len(atoms)*3,refcheck=False)
      
      self.pot = pot        # expects a depend object holding the potential
      self.f = f            # expects a depend object holding the force
      self.vir = vir        # expects a depend object holding the virial
      
      # use a compound object, since typically FF computes the force, potential and virial at once
      self.ufv = depend(name='ufv', func=self.get_all, value=self._ufv)
      for atom in self.atoms: 
         atom.q.add_dependant(self.ufv)
      
      self.ufv.add_dependant(self.pot); self.ufv.add_dependant(self.f); self.ufv.add_dependant(self.vir)
      self.pot._depend__func=self.get_pot; 
      self.f._depend__func=self.get_f; 
      self.vir._depend__func=self.get_vir;
      self.cell.V.add_dependant(self.vir)
      
   def get_all(self):
      """Dummy routine where no calculation is done"""

      return ( 0.0, numpy.zeros(len(self.atoms)*3), numpy.zeros((3,3),float) )
   
   def get_pot(self):
      """Calls get_all routine of forcefield to update potential"""

      (self._pot, self._f, self._vir) = self.ufv.get()
      return self._pot
     
   def get_f(self):
      """Calls get_all routine of forcefield to update force"""

      (self._pot, self._f, self._vir) = self.ufv.get()
      return self._f
   
   def get_vir(self):
      """Calls get_all routine of forcefield to update virial"""

      (self._pot, self._f, self._vir) = self.ufv.get()
      return self._vir/self.cell.V.get()
      

class pipeforce(forcefield):
   """Creates an interface between the force, potential and virial calculation
      and the wrapper. This uses an external program to do the update.
      Contains: _pot = potential, _f = forces, _vir = virial, _ufv = 
      (potential, forces, virial)
      Initialised by: ffield = pipeforce(dict)
      dict = {\"pipein\" : in, \"pipeout\" : out}
      in is a string containing the force input fifo, default = \"pipeforce\"
      out is a string containing the position output fifo, 
      default = \"pipepos\" """

   def __init__(self, pars=dict(pipein = "pipeforce", pipeout = "pipepos") ):
      super(pipeforce,self).__init__()
      self.pars = pars

      
   def bind(self, cell, atoms, pot, f, vir):
      """Binds the appropriate system objects to the force calculation, such
         that the potential, force and virial calculation automatically updates
         the same quantities in the system"""

      super(pipeforce,self).bind(cell=cell, atoms=atoms, pot=pot, f=f, vir=vir)
      self.ufv.func=self.get_all

   def get_all(self):
      """Updates _ufv using an external force calculator, using 
         named pipes/fifos"""

      self.fout=open(self.pars["pipeout"],"w")
      io_system.xml_write(self, self.fout)
      self.fout.close()
      self.fin=open(self.pars["pipein"],"r")
      force_vars = io_system.xml_read(self.fin)
      self.fin.close()
      return force_vars

class LJ(forcefield):
   """Creates an interface between the force, potential and virial calculation
      and the wrapper. This uses an internal Lennard-Jones potential
      to do the update.
      Contains: _pot = potential, _f = forces, _vir = virial, _ufv = 
      (potential, forces, virial)
      Initialised by: ffield = pipeforce(dict)
      dict = {\"eps\" : eps_value, \"sigma\" : sigma_value, \"rc\" : cutoff}
      eps_value = energy scale parameter, default = 1.0
      sigma value = length scale parameter, default = 1.0
      cutoff = cutoff radius, default = 2.5"""

   def __init__(self, pars=dict(eps = 1.0, sigma = 1.0, rc = 2.5) ):
      super(LJ,self).__init__() 
      self.eps = pars['eps']; self.sigma = pars['sigma'];  self.rc = pars['rc']
      self.vrc=4*self.eps*((self.sigma/self.rc)**12 - (self.sigma/self.rc)**6)
      
   def bind(self, cell, atoms, pot, f, vir):
      """Binds the appropriate system objects to the force calculation, such
         that the potential, force and virial calculation automatically updates
         the same quantities in the system"""

      super(LJ,self).bind(cell=cell, atoms=atoms, pot=pot, f=f, vir=vir)
      self.ufv.func=self.get_all

   def separation(self, atom_i, atom_j):
      """Calculates the vector and scalar separation between two atoms"""

      rij = self.cell.minimum_distance(atom_i, atom_j)
      r = math.sqrt(numpy.dot(rij, rij))
      return r, rij

   def LJ_fdf(self, r):
      """Calculates the force and potential at a given separation"""

      if (r > self.rc):
         return 0.0, 0.0
      else:
         sonr=self.sigma/r; sonr6=sonr**6
         return 4*self.eps*(6/r*sonr6*(2*sonr6-1)), 4*self.eps*(sonr6*(sonr6-1)) - self.vrc
         
   def LJ_fv(self, atom_i, atom_j):
      """Calculates the LJ potential and force between two atoms"""

      r, rij = self.separation(atom_i, atom_j)
      fij, v = self.LJ_fdf(r)
      fij/=r
      fij=rij*fij
      
      return fij, rij, v
      
   def get_all(self):
      """Updates _ufv using an internal LJ potential and force calculator"""

      natoms = len(self.atoms)

      _vir = numpy.zeros((3,3),float)
      _pot = 0.0
      _f = numpy.zeros(3*natoms,float)
   
      for i in range(natoms-1):
         atom_i = self.atoms[i]

         for j in range(i+1, natoms):
            atom_j = self.atoms[j]

            fij, rij, v = self.LJ_fv(atom_i, atom_j)
            _f[3*i:3*(i+1)]+=fij
            _f[3*j:3*(j+1)]-=fij
            _pot += v
            for k in range(3):
               for l in range(k, 3):
                  _vir[k,l] += fij[k]*rij[l]
            
      return (_pot, _f, _vir)
            
