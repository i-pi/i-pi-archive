import numpy as np
import math
from utils.depend import *
from utils.io import io_system
from driver.interface import Interface, RestartInterface
from utils.restart import *

class RestartForce(Restart):
   attribs = { "type" : (RestartValue,(str,"socket")) }
   fields =  { "interface" : (RestartInterface,()), "parameters" : (RestartValue, (dict,None)) }
   
   def store(self, force):
      if (type(force) is FFSocket):  
         self.type.store("socket")
         self.interface.store(force.socket)
         self.parameters.store(force.pars)
      else: self.type.store("unknown")
         

   def fetch(self):
      if self.type.fetch().upper() == "SOCKET": 
         force=FFSocket(pars=self.parameters.fetch(), interface=self.interface.fetch())
      else : force=ForceField()
      return force
      
class ForceField(dobject):
   """Creates an interface between the force, potential and virial calculation
      and the wrapper.
      Contains: _pot = potential, _f = forces, _vir = virial, _ufv = 
      [potential, forces, virial]
      Initialised by ffield = forcefield()"""

   def __init__(self):
      dset(self,"ufv", depend_value(name="ufv", func=self.get_all))
      
   def copy(self):    # creates a deep copy with everything but the bound bits 
      return type(self)()
      
   def bind(self, atoms, cell):
      self.atoms = atoms
      self.cell = cell
      dget(self,"ufv").add_dependency(dget(self.atoms,"q"))
      dget(self,"ufv").add_dependency(dget(self.cell,"h"))      
      dset(self,"pot",depend_value(name="pot", func=self.get_pot, dependencies=[dget(self,"ufv")] )  )
      
      fbase=np.zeros(atoms.natoms*3, float)
      dset(self,"f", depend_array(name="f", value=fbase, func=self.get_f, dependencies=[dget(self,"ufv")]) )
      # a bit messy, but we don't want to trigger quite yet the force calculation routine 
      dset(self,"fx", depend_array(name="fx", value=fbase[0:3*atoms.natoms:3]));
      dset(self,"fy", depend_array(name="fy", value=fbase[1:3*atoms.natoms:3]));
      dset(self,"fz", depend_array(name="fz", value=fbase[2:3*atoms.natoms:3]));
      depcopy(self,"f", self,"fx");      depcopy(self,"f", self,"fy");      depcopy(self,"f", self,"fz");
      dset(self,"vir", depend_array(name="vir", value=np.zeros((3,3),float),func=self.get_vir, dependencies=[dget(self,"ufv")] ) )

   def queue(self): pass
   
   def get_all(self):
      """Dummy routine where no calculation is done"""
      return [0.0, numpy.zeros(3*self.atoms.natoms), numpy.zeros((3,3),float)]

   def get_pot(self):
      """Calls get_all routine of forcefield to update potential"""

      [pot, f, vir] = self.ufv
      return pot

   def get_f(self):
      """Calls get_all routine of forcefield to update force"""

      [pot, f, vir] = self.ufv
      return f

   def get_vir(self):
      """Calls get_all routine of forcefield to update virial"""

      [pot, f, vir] = self.ufv
      vir[1,0]=0.0; vir[2,0:2]=0.0;
      return vir
      

import threading
class ForceBeads(dobject):
   """A class to collect many FF instances and parallelise getting the forces in a PIMD environment"""
   def __init__(self, beads=None, cell=None, force=None):
      if not (beads is None or cell is None or force is None): self.bind(beads, cell, force)
   
   def bind(self, beads, cell, force):
      self.natoms=beads.natoms
      self.nbeads=beads.nbeads

      self._forces=[];
      for b in range(self.nbeads):
         newf=force.copy()
         newf.bind(beads[b], cell)
         newf.blocking=False
         self._forces.append(newf)
               
      u=dget(self._forces[0],"f")
      dset(self,"f",depend_array(name="f",value=np.zeros((self.nbeads,3*self.natoms), float), func=self.f_gather,     
          dependencies=[dget(self._forces[b],"f")  for b in range(self.nbeads)] ) )
      dset(self,"pots",depend_array(name="pots", value=np.zeros(self.nbeads,float), func=self.pot_gather,     
          dependencies=[dget(self._forces[b],"pot")  for b in range(self.nbeads)] ) )
      dset(self,"pot",depend_value(name="pot", func=self.pot,     
          dependencies=[dget(self,"pots")] ) )
      dset(self,"fnm",depend_array(name="fnm",value=np.zeros((self.nbeads,3*self.natoms), float), func=self.b2nm_f, dependencies=[dget(self,"f")] ) )
      self.Cb2nm=beads.Cb2nm
      
   def b2nm_f(self): return np.dot(self.Cb2nm,depstrip(self.f))

   def _getbead(self, b, newf):
      newf[b]=self._forces[b].f
      return

   def queue(self):   
      for b in range(self.nbeads): self._forces[b].queue()

   def pot_gather(self): 
      self.queue()
      return np.array([b.pot for b in self._forces], float)
   def pot(self): return self.pots.sum()
      
   def f_gather(self): 
      start=time.time()
      newf=np.zeros((self.nbeads,3*self.natoms),float)
      
      self.queue()
#      print "time queueing", time.time()-start
#      print "update", self._forces[0].socket.time_update, 
#      print "distribute", self._forces[0].socket.time_distribute
#      self._forces[0].socket.time_update=0
#      self._forces[0].socket.time_distribute=0


      #serial
      for b in range(self.nbeads): newf[b]=self._forces[b].f
      # threaded      
#      bthreads=[]
#      print "starting threads"
#      for b in range(self.nbeads): 
#         thread=threading.Thread(target=self._getbead, args=(b,newf,))
#         thread.start()
#         bthreads.append(thread)

#      print "waiting threads"      
#      for b in range(self.nbeads): bthreads[b].join()
#      print "threads joined in"

      return newf
      
import time
class FFSocket(ForceField):

   def __init__(self, pars={}, interface=None):
      super(FFSocket,self).__init__() 
      if interface is None:
         self.socket=Interface()
      else:
         self.socket=interface
      self.pars=pars
      
      self.timer=0.0
      self.twall=0.0
      self.ncall=0
      self.request=None
   def copy(self):    # creates a deep copy with everything but the bound bits 
      return type(self)(self.pars, self.socket)

   def get_all(self):
      #print "computing forces"
      self.timer-= time.clock()
      self.twall-=time.time()
      
      if self.request is None: self.request=self.socket.queue(self.atoms, self.cell, self.pars)
      while self.request["status"] != "Done": time.sleep(self.socket.latency)
      self.socket.release(self.request)
      result=self.request["result"]
      self.request=None
      
      self.ncall+=1
      self.timer+=time.clock()
      self.twall+=time.time()      
      return result
      
   def queue(self):   
      if self.request is None and dget(self,"ufv").tainted():
         self.request=self.socket.queue(self.atoms, self.cell, self.pars)
      
      
#      self.timer-= time.clock()
#      self.timewait-= time.clock()
#      self.fout=open(self.pars["pipeout"],"w")      
#      self.timewait+=time.clock()
#      self.timewrite-=time.clock()
#      io_system.xml_write(self.atoms, self.cell, self.fout)      
#      self.timewrite+=time.clock()

#      self.fout.close()

#      self.timewait-= time.clock()
#      self.fin=open(self.pars["pipein"],"r")
#      self.timewait+= time.clock()
#      self.timeread-=time.clock()      
#      [pot, f, vir]=io_system.xml_read(self.fin)
#      self.timeread+=time.clock()
#      
#      self.fin.close()
#      self.timer+=time.clock()
#      self.ncall+=1      
#      return [pot, f, vir]

class FFLennardJones(ForceField):
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
      super(FFLennardJones,self).__init__() 
      self.eps = pars['eps']; self.sigma = pars['sigma'];  self.rc = pars['rc']
      self.vrc=4*self.eps*((self.sigma/self.rc)**12 - (self.sigma/self.rc)**6)

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

      print "Computing forces"
      natoms = self.atoms.natoms

      vir = numpy.zeros((3,3),float)
      pot = 0.0
      f = numpy.zeros(3*natoms,float)
   
      for i in range(natoms-1):
         atom_i = self.atoms[i]

         for j in range(i+1, natoms):
            atom_j = self.atoms[j]
            
            fij, rij, v = self.LJ_fv(atom_i, atom_j)
            f[3*i:3*(i+1)]+=fij
            f[3*j:3*(j+1)]-=fij
            pot += v
            for k in range(3):
               for l in range(k, 3):
                  vir[k,l] += fij[k]*rij[l]
            
      return [pot, f, vir]
