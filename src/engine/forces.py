import math, numpy
from utils.depend import *

class forcefield(object):
   _pot=0.0
   _f=numpy.zeros(0,float)
   _vir=numpy.zeros((3,3),float)
   
   def __init__(self, pars=dict()):
      self._ufv=(self._pot, self._f, self._vir)

   def bind(self, cell, atoms, pot, f, vir):
      self.cell=cell
      self.atoms=atoms      # expects a shallow copy of the system's atoms list            
      self._f.resize(len(atoms)*3,refcheck=False)
      
      self.pot = pot
      self.f = f                  # expecs a depend object holding the force
      self.vir = vir              # expecs a depend object holding the virial
      
      # use a compound object, since typically FF compute force, potential and virial at once
      self.ufv = depend(name='ufv', func=self.get_all, value=self._ufv)
      for at in self.atoms: at.q.add_dependant(self.ufv)
      
      self.ufv.add_dependant(self.pot); self.ufv.add_dependant(self.f); self.ufv.add_dependant(self.vir)
      self.pot._depend__func=self.get_pot; 
      self.f._depend__func=self.get_f; 
      self.vir._depend__func=self.get_vir;
      self.cell.V.add_dependant(self.vir)
      
   def get_all(self):
      return ( 0.0, numpy.zeros(len(self.atoms)*3), numpy.zeros((3,3),float) )
   
   def get_pot(self):
      (self._pot, self._f, self._vir) = self.ufv.get()
      return self._pot
     
   def get_f(self):
      (self._pot, self._f, self._vir) = self.ufv.get()
      return self._f
   
   def get_vir(self):
      (self._pot, self._f, self._vir) = self.ufv.get()
      return self._vir/self.cell.V.get()
      
      

class LJ(forcefield):
   def __init__(self, pars=dict(eps = 1.0, sigma = 1.0, rc = 2.5) ):
      super(LJ,self).__init__() 
      self.eps = pars['eps']; self.sigma = pars['sigma'];  self.rc = pars['rc']
      self.vrc=4*self.eps*((self.sigma/self.rc)**12 - (self.sigma/self.rc)**6)
      
   def bind(self, cell, atoms, pot, f, vir):
      super(LJ,self).bind(cell=cell, atoms=atoms, pot=pot, f=f, vir=vir)
      self.cell.V.add_dependant(self.ufv)
      self.ufv.func=self.get_all
      

   def separation(self, atom_i, atom_j):
      rij = self.cell.minimum_distance(atom_i, atom_j)
      r = math.sqrt(numpy.dot(rij, rij))
      return r, rij

   def LJ_force(self, r):
      if (r > self.rc):
         return 0.0
      else:
         return 4*self.eps*(12/r*(self.sigma/r)**12 - 6/r*(self.sigma/r)**6)

   def LJ_fij(self, atom_i, atom_j):
      fij = numpy.zeros(3, float)
      fji = numpy.zeros(3, float)
      r, rij = self.separation(atom_i, atom_j)
      f_tot = self.LJ_force(r)
      for i in range(3):
         fij[i] = f_tot*rij[i]/r
      fji = -fij
      return fij, fji, r, rij

   def LJ_pot(self, r):
      if (r > self.rc):
         return 0.0
      else:
         return 4*self.eps*((self.sigma/r)**12 - (self.sigma/r)**6 - (self.sigma/self.rc)**12 + (self.sigma/self.rc)**6)

   # a -- slightly -- faster LJ evaluator
   def LJ_fdf(self, r):
      if (r > self.rc):
         return 0.0, 0.0
      else:
         sonr=self.sigma/r; sonr6=sonr**6
         return 4*self.eps*(6/r*sonr6*(2*sonr6-1)), 4*self.eps*(sonr6*(sonr6-1)) - self.vrc
         
   def LJ_fv(self, atom_i, atom_j):
      r, rij = self.separation(atom_i, atom_j)
      fij, v = self.LJ_fdf(r)
      fij/=r
      fij=rij*fij
      
      return fij, rij, v
      
   def get_all(self):
#      print "computing all forces"
      natoms = len(self.atoms)

      _vir = numpy.zeros((3,3),float)
      _pot = 0.0
      _f = numpy.zeros(3*natoms,float)
   
      for i in range(natoms-1):
         atom_i = self.atoms[i]

         for j in range(i+1, natoms):
            atom_j = self.atoms[j]

#            fij, fji, r, rij = self.LJ_fij(atom_i, atom_j)
#            _f[3*i:3*(i+1)]+=fij
#            _f[3*j:3*(j+1)]+=fji
#            _pot += self.LJ_pot(r)
            fij, rij, v = self.LJ_fv(atom_i, atom_j)
            _f[3*i:3*(i+1)]+=fij
            _f[3*j:3*(j+1)]-=fij
            _pot += v
            for k in range(3):
               for l in range(k, 3):
                  _vir[k,l] += fij[k]*rij[l]
            
      return (_pot, _f, _vir)
            
