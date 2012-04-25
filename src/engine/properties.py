"""Holds the class which computes important properties of the system, and 
prepares them for output.

Classes:
   Properties: This is the class that holds all the algorithms to calculate
      the important properties that should be output.
   Trajectories: This class deals with outputting all position data in the
      appropriate format.
"""

__all__ = ['Properties', 'Trajectories']

import numpy as np
import math, random
from utils.depend import *
from utils.units import Constants
from utils.mathtools import h2abc
from utils.io import *
from atoms import *
from cell import *
from ensembles import *
from forces import *

_DEFAULT_FINDIFF = 1e-5
_DEFAULT_FDERROR = 1e-9
_DEFAULT_MINFID = 1e-12
class Properties(dobject):
   """A proxy to compute and output properties of the system.

   Takes the fundamental properties calculated during the simulation, and 
   prepares them for output. It also contains simple algorithms to calculate
   other properties not calculated during the simulation itself, so that 
   these can also be output.

   Attributes:
      fd_delta: A float giving the size of the finite difference
         parameter used in the Yamamoto kinetic energy estimator. Defaults 
         to _DEFAULT_FINDIFF.
      _DEFAULT_FDERROR: A float giving the size of the minimum precision 
         allowed for the finite difference calculation in the Yamamoto kinetic
         energy estimator.
      _DEFAULT_MINFID: A float giving the maximum displacement in the Yamamoto 
         kinetic energy estimator.
      dbeads: A dummy Beads object used in the Yamamoto kinetic energy
         estimator.
      dforces: A dummy Forces object used in the Yamamoto kinetic energy
         estimator.
      simul: The Simulation object containing the data to be output.
      ensemble: An ensemble object giving the objects necessary for producing
         the correct ensemble.
      beads: A beads object giving the atoms positions.
      cell: A cell object giving the system box.
      forces: A forcefield object giving the force calculator for each 
         replica of the system.
      property_dict: A dictionary containing all the properties that can be
         output.
      time: A float giving the time passed in the simulation.
      econs: A float giving the conserved quantity.
      kin: A float giving the classical kinetic energy estimator.
      pot: A float giving the potential energy estimator.
      temp: A float giving the classical kinetic temperature estimator.
      h: The elements of the unit cell matrix [h(x=vector_index,v=coordinate_index)]
      stress: An array giving the components of the classical stress tensor
         estimator.
      press: A float giving the classical pressure estimator.
      kin_cv: A float giving the quantum centroid virial kinetic energy 
         estimator.
      kstress_cv: An array giving the components of the quantum centroid virial
         kinetic stress tensor estimator.
      stress_cv: An array giving the components of the quantum centroid virial
         stress tensor estimator.
      press_cv: A float giving the quantum centroid virial pressure estimator.
      kin_yama: A float giving the quantum scaled coordinate estimator for the
         kinetic energy.
   """

   def __init__(self):
      """Initialises Properties."""

      self.property_dict = {}
      self.fd_delta = -_DEFAULT_FINDIFF
      
   def bind(self, simul):
      """Binds the necessary objects from the simulation to calculate the
      required properties.

      This function takes the appropriate simulation object, and creates the
      property_dict object which holds all the objects which can be output.
      It is given by: 
      {'time': Time elapsed,
      'step': The current time step,
      'conserved': Conserved quantity,
      'temperature': Classical kinetic temperature estimator,
      'volume': Simulation box volume,
      'h': Cell vector matrix. Requires arguments x and v to give h[x,v],
      'potential': Potential energy estimator,
      'kinetic_md': Classical kinetic energy estimator,
      'kinetic_cv': Quantum centroid virial kinetic energy estimator,
      'stress_md': The classical stress tensor estimator. Requires arguments
         x and v, to give stress[x,v],
      'pressure_md': Classical pressure estimator,
      'stress_cv': The quantum centroid virial estimator of 
         the stress tensor. Requires arguments x and v, to give stress[x,v],
      'pressure_cv': Quantum centroid virial pressure estimator,
      'kstress_cv': Quantum centroid virial kinetic stress tensor estimator.
         Requires arguments x and v, to give kstress[x,v],
      'kin_yama': Quantum scaled coordinate kinetic energy estimator,
      'linlin': The scaled Fourier transform of the momentum distribution.
         Given by n(x) in Lin Lin et al., Phys. Rev. Lett. 105, 110602.}.

      Args:
         simul: The Simulation object to be bound.
      """

      self.ensemble = simul.ensemble
      self.beads = simul.beads
      self.cell = simul.cell
      self.forces = simul.forces
      self.simul = simul


      self.property_dict["step"] = lambda: (1 + self.simul.step)
      self.property_dict["time"] = lambda: (1 + self.simul.step)*self.ensemble.dt
      self.property_dict["conserved"] = self.get_econs
      self.property_dict["temperature"] = self.get_temp
      self.property_dict["volume"] = lambda: self.cell.V
      
      self.property_dict["h"] = self.wrap_cell
      
      self.property_dict["potential"] = lambda: self.forces.pot/self.beads.nbeads
      self.property_dict["kinetic_md"] = lambda: self.beads.kin/self.beads.nbeads
      self.property_dict["kinetic_cv"] = self.get_kincv

      self.property_dict["stress_md"] = self.get_stress
      self.property_dict["pressure_md"] = self.get_press
      self.property_dict["stress_cv"] = self.get_stresscv
      self.property_dict["pressure_cv"] = self.get_presscv
      self.property_dict["kstress_cv"] = self.get_kstresscv

      self.property_dict["kin_yama"] = self.get_kinyama

      self.property_dict["linlin"] = self.get_linlin

      self.dbeads = simul.beads.copy()
      self.dforces = ForceBeads()
      self.dforces.bind(self.dbeads, self.simul.cell,  self.simul._forcemodel)
      
   def __getitem__(self, key):
      """Retrieves the item given by key.

      Note that if the key contains a string (arg1=value1; arg2=value2; ... )
      then it will add the appropriate arguments and value pairs
      to the calculation function of the property.

      Args:
         key: A string contained in property_dict.

      Returns:
         The property labelled by the keyword key.
      """

      args = {}
      if '(' in key:
         # If the property has additional arguments
         argstart = key.find('(')
         argstop = key.find(')', argstart)
         if argstop == -1:
            raise ValueError("Incorrect format in property name " + key)
         
         argstr = key[argstart:argstop+1]
         key = key[0:argstart] # strips the arguments from key name

         arglist = io_xml.read_dict(argstr, delims="()", split=";", key_split="=")
         
         return self.property_dict[key](**arglist)
      else:
         return self.property_dict[key]()

   def get_temp(self):
      """Calculates the classical kinetic temperature estimator.

      Note that in the case that the centre of mass constraint there will be
      one less degree of freedom than without, so this has to be taken into
      account when calculating the kinetic temperature.
      """

      if self.ensemble.fixcom:
         mdof=3 
      else:
         mdof=0
      return self.beads.kin/(0.5*Constants.kb*(3*self.beads.natoms*self.beads.nbeads - mdof)*self.beads.nbeads)

   def get_econs(self):
      """Calculates the conserved quantity estimator."""

      return self.ensemble.econs/(self.beads.nbeads*self.beads.natoms)

   def get_stress(self, x=0, v=0):
      """Calculates the classical kinetic energy estimator.

      Returns stress[x,v].
      """

      x = int(x)
      v = int(v)
      stress = (self.forces.vir + self.beads.kstress)/self.cell.V
      return stress[x,v]

   def get_press(self):
      """Calculates the classical pressure estimator."""

      stress = (self.forces.vir + self.beads.kstress)/self.cell.V
      return np.trace(stress)/3.0

   def kstress_cv(self):
      """Calculates the quantum central virial kinetic stress tensor 
      estimator.
      """

      kst = np.zeros((3,3),float)
      q = depstrip(self.beads.q)
      qc = depstrip(self.beads.qc)
      na3 = 3*self.beads.natoms;
      for b in range(self.beads.nbeads):
         for i in range(3):
            for j in range(i,3):
               kst[i,j] += np.dot(q[b,i:na3:3] - qc[i:na3:3], depstrip(self.forces.f[b])[j:na3:3])

      kst *= -1/self.beads.nbeads
      for i in range(3):
         kst[i,i] += Constants.kb*self.ensemble.temp*(3*self.beads.natoms) 
      return kst

   def get_kstresscv(self, x=0, v=0):        
      """Calculates the quantum central virial kinetic stress tensor 
      estimator.

      Returns kstress[x,v].
      """

      x = int(x)
      v = int(v)
      return self.kstress()[x,v]

   def get_stresscv(self, x=0, v=0):
      """Calculates the quantum central virial stress tensor estimator.

      Returns stress[x,v].
      """

      x = int(x)
      v = int(v)
      stress = (self.forces.vir + self.kstress())/self.cell.V                  
      return stress[x,v]

   def get_presscv(self):
      """Calculates the quantum central virial pressure estimator."""

      return np.trace(self.forces.vir + self.kstress())/(3.0*self.cell.V)
   
   def get_kincv(self):        
      """Calculates the quantum central virial kinetic energy estimator."""

      kcv=0.0
      for b in range(self.beads.nbeads):
         kcv += np.dot(depstrip(self.beads.q[b]) - depstrip(self.beads.qc), depstrip(self.forces.f[b]))
      kcv *= -0.5/self.beads.nbeads
      kcv += 0.5*Constants.kb*self.ensemble.temp*(3*self.beads.natoms) 
      return kcv

   def get_kinyama(self):              
      """Calculates the quantum scaled coordinate kinetic energy estimator.

      Uses a finite difference method to calculate the kinetic energy estimator
      without requiring the forces as for the centroid virial estimator.
      """
      
      dbeta = abs(self.fd_delta)
      
      v0 = self.forces.pot/self.beads.nbeads
      while True: 
         splus = math.sqrt(1.0 + dbeta)
         sminus = math.sqrt(1.0 - dbeta)
         
         for b in range(self.beads.nbeads):
            self.dbeads[b].q = self.beads.centroid.q*(1.0 - splus) + splus*self.beads[b].q
         vplus = self.dforces.pot/self.beads.nbeads
         
         for b in range(self.beads.nbeads):
            self.dbeads[b].q = self.beads.centroid.q*(1.0 - sminus) + sminus*self.beads[b].q      
         vminus = self.dforces.pot/self.beads.nbeads

         kyama = ((1.0 + dbeta)*vplus - (1.0 - dbeta)*vminus)/(2*dbeta) - v0
         kyama += 0.5*Constants.kb*self.ensemble.temp*(3*self.beads.natoms) 
         if (self.fd_delta < 0 and abs((vplus + vminus)/(v0*2) - 1.0) > _DEFAULT_FDERROR and dbeta > _DEFAULT_MINFID):
            dbeta *= 0.5
            print "Reducing displacement in Yamamoto kinetic estimator"
            continue
         else:
            break
         
      return kyama

   def opening(self, bead):
      """Path opening function.

      Used in the Lin Lin momentum distribution estimator.
      """

      return bead/self.beads.nbeads + 0.5*(1.0/self.beads.nbeads - 1.0)

   def get_linlin(self, ux=0, uy=0, uz=0, atom=0):
      """Gives the estimator for the momentum distribution, by opening the 
      ring polymer path.

      Args:
         ux: The x component of the opening vector.
         uy: The y component of the opening vector.
         uz: The z component of the opening vector.
         atom: The atom for which the path will be opened.
      """

      u = np.array([float(ux), float(uy), float(uz)])
      atom = int(atom)
      for at in range(self.beads.natoms):
         if at == atom:
            for bead in range(self.beads.nbeads):
               self.dbeads.q[bead,3*at:3*(at+1)] = self.opening(bead)*u + self.beads.q[bead,3*at:3*(at+1)]
         else:
            self.dbeads.q[:,3*at:3*(at+1)] = self.beads.q[:,3*at:3*(at+1)]

      n0 = math.exp(self.beads.m[atom]*np.dot(u,u)*self.ensemble.temp*Constants.kb/(2*Constants.hbar**2))

      return n0*math.exp(-(self.dforces.pot - self.forces.pot)/(self.ensemble.ntemp*Constants.kb))

   def wrap_cell(self, x=0, v=0):
      """Returns the the x-th component of the v-th cell vector."""
   
      return self.cell.h[x,v]


class Trajectories(dobject):
   """A simple class to take care of output of trajectory data.

   Attributes:
      format: The file format for the output files.
      simul: The simulation object from which the position data will be 
         obtained.
      fatom: A dummy beads object used so that individual replica trajectories
         can be output.
   """
   
   def __init__(self, format = "pdb"):
      """Initialises a Trajectories object.

      Args:
         format: A format string giving the file format for the output files.
            Defaults to 'pdb'.
      """

      self.format = format
      
   def bind(self, simul):
      """ Binds to a simulation object to fetch atomic and force data.
      
      Args:
         simul: The simulation object that will be managed by this Trajectories.
      """

      self.simul = simul
      self.fatom = simul.beads[0].copy()
      
      # a few, "fancier", per-atom properties
      dset(self, "atomic_kincv", depend_array(name="atomic_kincv", value=np.zeros(self.simul.beads.natoms*3),
           func=self.get_akcv, dependencies=[dget(self.simul.forces,"f"), dget(self.simul.beads,"q"), dget(self.simul.beads,"qc"), dget(self.simul.ensemble,"temp")]))      
      dset(self, "atomic_kod", depend_array(name="atomic_kod", value=np.zeros(self.simul.beads.natoms*3),
           func=self.get_akcv_od, dependencies=[dget(self.simul.forces,"f"), dget(self.simul.beads,"q"), dget(self.simul.beads,"qc"), dget(self.simul.ensemble,"temp")]))      
      
      
   def get_akcv(self):
      """Calculates the contribution to the kinetic energy due to each degree
      of freedom.
      """

      rv = np.zeros(self.simul.beads.natoms*3)
      for b in range(self.simul.beads.nbeads):
         rv[:] += (self.simul.beads.q[b]-self.simul.beads.qc)*self.simul.forces.f[b]
      rv *= -0.5/self.simul.nbeads
      rv += 0.5*Constants.kb*self.simul.ensemble.temp
      return rv

   def get_akcv_od(self):
      """Calculates the "off-diagonal" contribution to the kinetic energy tensor 
      due to each atom.
      """

      rv = np.zeros((self.simul.beads.natoms,3))
      # helper arrays to make it more transparent what we are computing
      dq = np.zeros((self.simul.beads.natoms,3))
      f = np.zeros((self.simul.beads.natoms,3))
      for b in range(self.simul.beads.nbeads):
         dq[:] = (self.simul.beads.q[b]-self.simul.beads.qc).reshape((self.simul.beads.natoms,3))
         f[:] = self.simul.forces.f[b].reshape((self.simul.beads.natoms,3))
         rv[:,0] += dq[:,0]*f[:,1]+dq[:,1]*f[:,0]
         rv[:,1] += dq[:,1]*f[:,2]+dq[:,2]*f[:,1]
         rv[:,2] += dq[:,0]*f[:,2]+dq[:,2]*f[:,0]
      rv *= 0.5
      rv *= -0.5/self.simul.nbeads
      # rv += 0.5*Constants.kb*self.simul.ensemble.temp
      
      return rv.reshape(self.simul.beads.natoms*3)
         
   def print_traj(self, what, stream, b=0):
      """Prints out a frame of a trajectory for the specified quantity and bead.

      Args:
         what: A string specifying what to print.
         b: The bead index. Defaults to 0.
         stream: A reference to the stream on which data will be printed.
      """
      
      if what == "positions":
         self.fatom.q = self.simul.beads.q[b]
      elif what == "velocities":
         self.fatom.q = self.simul.beads.p[b]/self.simul.beads.m3[0]      
      elif what == "forces":
         self.fatom.q = self.simul.forces.f[b]
      elif what == "kinetic_cv":
         self.fatom.q = self.atomic_kincv 
      elif what == "kodterms_cv":
         self.fatom.q = self.atomic_kod
      elif what == "centroid":
         self.fatom.q = self.simul.beads.qc
      else:
         raise IndexError("<" + what + "> is not a recognized trajectory output")
      
      if self.format == "pdb":
         io_pdb.print_pdb(self.fatom, self.simul.cell, stream, title=("Step:  %10d  Bead:   %5d " % (self.simul.step+1, b) ) )
      elif self.format == "xyz":
         io_xyz.print_xyz(self.fatom, self.simul.cell, stream, title=("Step:  %10d  Bead:   %5d " % (self.simul.step+1, b) ) )
