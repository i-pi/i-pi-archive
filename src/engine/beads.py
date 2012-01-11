import numpy as np
import math
from utils.depend import *
from utils.restart import *
from engine.atoms import Atoms

from utils import units
         
class Beads(dobject):
   """Storage for the atoms positions and velocities. Everything is stored as 3*n-sized contiguous arrays,
   and a convenience-access is provided through a list of Atom objects"""   

   def __init__(self, natoms, nbeads):
      self.natoms=natoms
      self.nbeads=nbeads

      # assumes masses and names are same on different beads
      dset(self,"names",depend_array(name="names",value=np.zeros(natoms, np.dtype('|S6'))) )            
      dset(self,"m",depend_array(name="m",value=np.zeros(natoms, float)) )     
      # Interface to get a 3*n-sized array with masses      
      dset(self,"m3",depend_array(name="m3",value=np.zeros((nbeads,3*natoms), float),func=self.mtom3, dependencies=[dget(self,"m")]))
      dset(self,"sm3",depend_array(name="sm3",value=np.zeros((nbeads,3*natoms), float),func=self.m3tosm3, dependencies=[dget(self,"m3")]))
            
      # sets up synchronized storage for normal-modes and beads representations
      sync_q=synchronizer();    sync_p=synchronizer()
      dset(self,"q",   depend_array(name="q",   value=np.zeros((nbeads,3*natoms), float), func={"qnm":self.nm2b_q}, synchro=sync_q ) )
      dset(self,"p",   depend_array(name="p",   value=np.zeros((nbeads,3*natoms), float), func={"pnm":self.nm2b_p}, synchro=sync_p ) )
      dset(self,"qnm", depend_array(name="qnm", value=np.zeros((nbeads,3*natoms), float), func={"q":self.b2nm_q},   synchro=sync_q ) )
      dset(self,"pnm", depend_array(name="pnm", value=np.zeros((nbeads,3*natoms), float), func={"p":self.b2nm_p},   synchro=sync_p ) )

      
      # sets up matrices for normal modes transformation
      self.Cb2nm=np.zeros((nbeads,nbeads))
      self.Cb2nm[0,:]=math.sqrt(1.0/nbeads)
      for i in range(1,nbeads/2):
         for j in range(nbeads):  self.Cb2nm[i,j]=math.sqrt(2.0/nbeads)* math.cos(2*math.pi*j*i/float(nbeads))
      if (nbeads%2)==0:  
         self.Cb2nm[nbeads/2,0:nbeads:2]=math.sqrt(1.0/nbeads); self.Cb2nm[nbeads/2,1:nbeads:2]=-math.sqrt(1.0/nbeads)
      for i in range(nbeads/2+1,nbeads):
         for j in range(nbeads):  self.Cb2nm[i,j]=math.sqrt(2.0/nbeads)* math.sin(2*math.pi*j*i/float(nbeads))

      self.Cnm2b= self.Cb2nm.T.copy()

      #  Access to individual beads' properties via Atoms objects     
      self._blist=[ Atoms(natoms, _prebind=( self.q[i,:], self.p[i,:], self.m,  self.names )) for i in range(nbeads) ]

      dset(self,"qc",depend_array(name="qc",value=np.zeros(3*natoms, float), func=self.get_qc, dependencies=[dget(self,"qnm")] ) )      
      dset(self,"pc",depend_array(name="pc",value=np.zeros(3*natoms, float), func=self.get_pc, dependencies=[dget(self,"pnm")] ) )      
      self.centroid = Atoms(natoms, _prebind=(self.qc, self.pc, self.m, self.names))
      
      dset(self,"vpath",depend_value(name="vpath", func=self.vpath, dependencies=[dget(self,"q")]) )
      dset(self,"fpath",depend_array(name="fpath", value=np.zeros((nbeads,3*natoms), float), func=self.fpath, dependencies=[dget(self,"q")]) )
      dset(self,"kins",depend_array(name="kins",value=np.zeros(nbeads, float), func=self.kin_gather, dependencies=[dget(b,"kin") for b in self._blist] ) )
      dset(self,"kin",depend_value(name="kin", func=self.get_kin, dependencies=[dget(self,"kins")]) )
      dset(self,"kstress",depend_array(name="kstress",value=np.zeros((3,3), float), func=self.get_kstress, dependencies=[dget(b,"kstress") for b in self._blist] ) )

   def copy(self):
      newbd=Beads(self.natoms, self.nbeads)
      newbd.q[:]=self.q
      newbd.p[:]=self.p
      newbd.m[:]=self.m
      newbd.names[:]=self.names
      return newbd

   def m3tosm3(self): return np.sqrt(self.m3)
   def mtom3(self): m3=np.zeros((self.nbeads,3*self.natoms),float); m3[:,0:3*self.natoms:3]=self.m; m3[:,1:3*self.natoms:3]=m3[:,0:3*self.natoms:3]; m3[:,2:3*self.natoms:3]=m3[:,0:3*self.natoms:3]; return m3
   
   def nm2b_q(self): return np.dot(self.Cnm2b,depstrip(self.qnm))
   def nm2b_p(self): return np.dot(self.Cnm2b,depstrip(self.pnm))
   def b2nm_q(self): return np.dot(self.Cb2nm,depstrip(self.q))
   def b2nm_p(self): return np.dot(self.Cb2nm,depstrip(self.p))
   def get_qc(self): return depstrip(self.qnm)[0,:]/math.sqrt(self.nbeads)
   def get_pc(self): return depstrip(self.pnm)[0,:]/math.sqrt(self.nbeads)
   
   def kin_gather(self):  return np.array([b.kin for b in self._blist])
   def get_kin(self):     return self.kins.sum()
   def get_kstress(self):     
      ks=np.zeros((3,3),float)
      for b in self: ks+=b.kstress
      return ks
      
   def vpath(self):  # this is actually the path harmonic potential without the multiplication by omega_n^2
      epath=0.0
      q=depstrip(self.q)
      m=depstrip(self.m3[0])
      for b in range(self.nbeads):
         if b>0 : dq=q[b,:]-q[b-1,:]
         else:    dq=q[b,:]-q[self.nbeads-1,:]
         epath+=np.dot(dq, m*dq)
      return epath*0.5  

   def fpath(self):  # this is actually the path harmonic force without the multiplication by omega_n^2
      nbeads=self.nbeads; natoms=self.natoms
      f=np.zeros((nbeads,3*natoms),float)
      
      q=depstrip(self.q)
      m=depstrip(self.m3[0])
      for b in range(nbeads):
         if b>0 : dq=q[b,:]-q[b-1,:]
         else:    dq=q[b,:]-q[self.nbeads-1,:]
         dq*=m
         f[b]-=dq
         if b>0: f[b-1]+=dq
         else:   f[nbeads-1]+=dq
      return f
      
   def __len__(self): return self.nbeads
   
   def __getitem__(self,index):
      return self._blist[index]

   def __setitem__(self,index,value):
      self._blist[index].p[:]=value.p
      self._blist[index].q[:]=value.q 
      self._blist[index].m[:]=value.m      
      self._blist[index].names[:]=value.names         
   
                        
from utils.io.io_pdb import read_pdb      
class RestartBeads(Restart):
   fields={ "nbeads" : (RestartValue, (int, 0)), "natoms" : (RestartValue, (int, 0)),  "q" : (RestartArray,(float,np.zeros(0))),  "p" : (RestartArray,(float,np.zeros(0))),
            "m" : (RestartArray,(float, np.zeros(0))),  "names" : (RestartArray,(str,np.zeros(0, np.dtype('|S6')))),
            "init_temp": (RestartValue, (float, -1.0))  }
   
   def __init__(self, beads=None):
      super(RestartBeads,self).__init__()
      if not beads is None: self.store(beads)
                       
   def store(self, beads):
      self.natoms.store(beads.natoms)
      self.nbeads.store(beads.nbeads)

      self.q.store(depstrip(beads.q))
      self.p.store(depstrip(beads.p))
      self.m.store(depstrip(beads.m))
      self.names.store(depstrip(beads.names))
      
   def fetch(self):
      beads=Beads(self.natoms.fetch(),self.nbeads.fetch())
      beads.q=self.q.fetch()
      beads.p=self.p.fetch()      
      beads.m=self.m.fetch()   
      beads.names=self.names.fetch()
      return beads
