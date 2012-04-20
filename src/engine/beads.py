"""Contains the classes which deal with all the beads.

Used for holding information about the beads, including their positions, masses
momenta and kinetic energy. Has different objects for the position and normal
mode representations, and has a special centroid atoms object for when the 
centroid coordinate is required.

Classes:
   Beads: Class with methods dealing with all the beads.
"""

__all__ = ['Beads']

import numpy as np
import math
from utils.depend import *
from engine.atoms import Atoms
from utils import units
         
class Beads(dobject):
   """Storage for the beads positions and velocities.

   Everything is stored as (nbeads,3*natoms) sized contiguous arrays,
   and a convenience-access to each replica of the system is provided through a
   list of Atoms objects. Contains arrays of both the normal mode representation
   and the position representation, and various sized arrays for the atom
   labels and masses. Also contains the potential and force between 
   neighbouring replicas.

   Attributes:
      natoms: The number of atoms.
      nbeads: The number of beads.
      Cb2nm: Transformation matrix between the bead and normal mode
         representations
      Cnm2b: Transformation matrix between the normal mode and bead
         representations
      _blist: A list of Atoms objects for each replica of the system. Each 
         replica is assumed to have the same mass and atom label.
      centroid: An atoms object giving the centroid coordinate of the beads.

   Depend objects:
      names: An array giving the atom names.
      m: An array giving the atom masses.
      m3: An array giving all the bead masses.
      sm3: An array giving the square root of all the bead masses.
      q: An array giving all the bead positions.
      p: An array giving all the bead momenta.
      qnm: An array giving the normal mode representation of the beads.
      pnm: An array giving the normal mode representation of the bead momenta.
      qc: An array giving the centroid positions.
      pc: An array giving the centroid momenta.
      vpath: The spring potential between the beads, divided by omegan**2.
      fpath: The spring force between the beads, divided by omegan**2.
      kins: A list of the kinetic energy of each replica.
      kin: The total kinetic energy of the system. Note that this is not the
         same as the estimate of the kinetic energy of the system, which is
         contained in the properties module. 
      kstress: The total kinetic stress tensor for the system.
   """   

   def __init__(self, natoms, nbeads):
      """Initialises Beads.

      Args:
         natoms: Number of atoms.
         nbeads: Number of beads.
      """

      self.natoms = natoms
      self.nbeads = nbeads

      dset(self,"names",depend_array(name="names",value=np.zeros(natoms, np.dtype('|S6'))) )            
      dset(self,"m",depend_array(name="m",value=np.zeros(natoms, float)) )     
      dset(self,"m3",depend_array(name="m3",value=np.zeros((nbeads,3*natoms), float),func=self.mtom3, dependencies=[dget(self,"m")]))
      dset(self,"sm3",depend_array(name="sm3",value=np.zeros((nbeads,3*natoms), float),func=self.m3tosm3, dependencies=[dget(self,"m3")]))
            
      sync_q = synchronizer()
      sync_p = synchronizer()
      dset(self,"q",depend_array(name="q",value=np.zeros((nbeads,3*natoms), float), func={"qnm":self.nm2b_q}, synchro=sync_q) )
      dset(self,"p",depend_array(name="p",value=np.zeros((nbeads,3*natoms), float), func={"pnm":self.nm2b_p}, synchro=sync_p) )
      dset(self,"qnm",depend_array(name="qnm",value=np.zeros((nbeads,3*natoms), float), func={"q":self.b2nm_q}, synchro=sync_q) )
      dset(self,"pnm",depend_array(name="pnm",value=np.zeros((nbeads,3*natoms), float), func={"p":self.b2nm_p}, synchro=sync_p) )
      
      self.Cb2nm = np.zeros((nbeads,nbeads))
      self.Cb2nm[0,:] = math.sqrt(1.0/nbeads)
      for i in range(1,nbeads/2+1):
         for j in range(nbeads):
            self.Cb2nm[i,j] = math.sqrt(2.0/nbeads)*math.cos(2*math.pi*j*i/float(nbeads))
      if (nbeads%2) == 0:  
         self.Cb2nm[nbeads/2,0:nbeads:2] = math.sqrt(1.0/nbeads)
         self.Cb2nm[nbeads/2,1:nbeads:2] = -math.sqrt(1.0/nbeads)
      for i in range(nbeads/2+1, nbeads):
         for j in range(nbeads):
            self.Cb2nm[i,j] = math.sqrt(2.0/nbeads)*math.sin(2*math.pi*j*i/float(nbeads))

      self.Cnm2b = self.Cb2nm.T.copy()

      self._blist = [Atoms(natoms, _prebind=( self.q[i,:], self.p[i,:], self.m,  self.names )) for i in range(nbeads) ]

      dset(self,"qc",depend_array(name="qc",value=np.zeros(3*natoms, float), func=self.get_qc, dependencies=[dget(self,"qnm")] ) )      
      dset(self,"pc",depend_array(name="pc",value=np.zeros(3*natoms, float), func=self.get_pc, dependencies=[dget(self,"pnm")] ) )      
      self.centroid = Atoms(natoms, _prebind=(self.qc, self.pc, self.m, self.names))
      
      dset(self,"vpath",depend_value(name="vpath", func=self.vpath, dependencies=[dget(self,"q")]) )
      dset(self,"fpath",depend_array(name="fpath", value=np.zeros((nbeads,3*natoms), float), func=self.fpath, dependencies=[dget(self,"q")]) )
      dset(self,"kins",depend_array(name="kins",value=np.zeros(nbeads, float), func=self.kin_gather, dependencies=[dget(b,"kin") for b in self._blist] ) )
      dset(self,"kin",depend_value(name="kin", func=self.get_kin, dependencies=[dget(self,"kins")]) )
      dset(self,"kstress",depend_array(name="kstress",value=np.zeros((3,3), float), func=self.get_kstress, dependencies=[dget(b,"kstress") for b in self._blist] ) )

   def copy(self):
      """Creates a new beads object from the original.

      Returns:
         A Beads object with the same q, p, m and names arrays as the original.
      """

      newbd = Beads(self.natoms, self.nbeads)
      newbd.q[:] = self.q
      newbd.p[:] = self.p
      newbd.m[:] = self.m
      newbd.names[:] = self.names
      return newbd

   def m3tosm3(self):
      """Takes the mass array and returns the square rooted mass array."""

      return np.sqrt(self.m3)

   def mtom3(self):
      """Takes the mass array for each bead and returns one with an element
      for each degree of freedom.

      Returns:
         An array of size (nbeads,3*natoms), with each element corresponding
         to the mass associated with the appropriate degree of freedom in q.
      """

      m3 = np.zeros((self.nbeads,3*self.natoms),float)
      m3[:,0:3*self.natoms:3] = self.m
      m3[:,1:3*self.natoms:3] = m3[:,0:3*self.natoms:3]
      m3[:,2:3*self.natoms:3] = m3[:,0:3*self.natoms:3]
      return m3
   
   def nm2b_q(self):
      """Takes the normal mode representation for q and returns the bead 
      representation.
      """

      return np.dot(self.Cnm2b,depstrip(self.qnm))

   def nm2b_p(self):
      """Takes the normal mode representation for p and returns the bead 
      representation.
      """

      return np.dot(self.Cnm2b,depstrip(self.pnm))

   def b2nm_q(self):
      """Takes the bead representation for q and returns the normal mode 
      representation.
      """

      return np.dot(self.Cb2nm,depstrip(self.q))

   def b2nm_p(self):
      """Takes the bead representation for p and returns the normal mode
      representation.
      """

      return np.dot(self.Cb2nm,depstrip(self.p))

   def get_qc(self):
      """Gets the centroid coordinates."""

      return depstrip(self.qnm)[0,:]/math.sqrt(self.nbeads)

   def get_pc(self):
      """Gets the centroid momenta."""

      return depstrip(self.pnm)[0,:]/math.sqrt(self.nbeads)
   
   def kin_gather(self):
      """Gets the kinetic energy for all the replicas.

      Returns:
         A list of the kinetic energy for each system.
      """

      return np.array([b.kin for b in self._blist])

   def get_kin(self):
      """Gets the total kinetic energy of all the replicas.

      Note that this does not correspond to the total kinetic energy estimate
      for the system.

      Returns:
         The sum of the kinetic energy of each replica.
      """

      return self.kins.sum()

   def get_kstress(self):     
      """Calculates the total kinetic stress tensor of all the replicas.

      Note that this does not correspond to the total kinetic stress tensor
      estimate for the system.

      Returns:
         The sum of the kinetic stress tensor of each replica.
      """

      ks = np.zeros((3,3),float)
      for b in self:
         ks += b.kstress
      return ks
      
   def vpath(self):
      """Calculates the spring potential between the replicas.

      Note that this is actually the harmonic potential without being
      multiplied by the factor omegan**2, which is only available in the
      ensemble as the temperature is required to calculate it.
      """

      epath = 0.0
      q = depstrip(self.q)
      m = depstrip(self.m3[0])
      for b in range(self.nbeads):
         if b > 0:
            dq = q[b,:] - q[b-1,:]
         else:
            dq = q[b,:] - q[self.nbeads-1,:]
         epath += np.dot(dq, m*dq)
      return epath*0.5  

   def fpath(self):
      """Calculates the spring force between the replicas.

      Note that this is actually the harmonic force without being
      multiplied by the factor omegan**2, which is only available in the
      ensemble as the temperature is required to calculate it.
      """

      nbeads = self.nbeads
      natoms = self.natoms
      f = np.zeros((nbeads,3*natoms),float)
      
      q = depstrip(self.q)
      m = depstrip(self.m3[0])
      for b in range(nbeads):
         if b > 0:
            dq = q[b,:] - q[b-1,:]
         else:
            dq = q[b,:] - q[self.nbeads-1,:]
         dq *= m
         f[b] -= dq
         if b > 0:
            f[b-1] += dq
         else:
            f[nbeads-1] += dq
      return f
      
   def __len__(self):
      """Length function.

      This is called whenever the standard function len(beads) is used.

      Returns:
         The number of beads.
      """

      return self.nbeads
   
   def __getitem__(self,index):
      """Overwrites standard getting function.

      This is called whenever the standard function beads[index] is used.
      Returns an Atoms object with the appropriate position and momenta arrays.

      Args:
         index: The index of the replica of the system to be accessed.

      Returns:
         The replica of the system given by the index.
      """

      return self._blist[index]

   def __setitem__(self,index,value):
      """Overwrites standard setting function.

      This is called whenever the standard function beads[index]=value is used.
      Changes the position and momenta of the appropriate slice of the global
      position and momentum arrays to those given by value.

      Args:
         index: The replica of the system to be changed.
         value: The Atoms object that holds the new values.
      """

      self._blist[index].p[:] = value.p
      self._blist[index].q[:] = value.q 
      self._blist[index].m[:] = value.m      
      self._blist[index].names[:] = value.names         
