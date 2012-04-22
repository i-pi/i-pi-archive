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
      cell_params: A list giving lattice vector lengths and the angles
         between them.
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
      'conserved': Conserved quantity,
      'kinetic_md': Classical kinetic energy estimator,
      'potential': Potential energy estimator,
      'temperature': Classical kinetic temperature estimator,
      'cell_parameters': Lattice vector lengths and the angles between them,
      'volume': Simulation box volume,
      'stress_md.xx': The xx component of the classical stress tensor estimator,
      'pressure_md': Classical pressure estimator
      'kinetic_cv': Quantum centroid virial kinetic energy estimator,
      'stress_cv.xx': xx component of the quantum centroid virial estimator of 
         the stress tensor,
      'pressure_cv': Quantum centroid virial pressure estimator
      'kinetic_yamamoto': Quantum scaled coordinate kinetic energy estimator}.

      Args:
         simul: The Simulation object to be bound.
      """

      self.ensemble = simul.ensemble
      self.beads = simul.beads
      self.cell = simul.cell
      self.forces = simul.forces
      self.simul = simul      

      self.add_property(prop_name="step", dep_name="step", func=self.get_step, dependencies=[dget(self.simul, "step")])
      self.add_property(prop_name="time", dep_name="time", func=self.get_time, dependencies=[dget(self.simul, "step"), dget(self.ensemble, "dt")])
      self.add_property(prop_name="conserved", dep_name="econs", func=self.get_econs, dependencies=[dget(self.ensemble, "econs")])
      self.add_property(prop_name="kinetic_md", dep_name="kin", func=self.get_kin, dependencies=[dget(self.beads, "kin"), dget(self.cell, "kin")])
      self.add_property(prop_name="potential", dep_name="pot", func=self.get_pot, dependencies=[dget(self.forces, "pot")])
      self.add_property(prop_name="temperature", dep_name="temp", func=self.get_temp, dependencies=[dget(self.beads, "kin")])

      self.property_dict["volume"] = dget(self.cell,"V")
      self.add_property(prop_name="cell_parameters", dep_name="cell_params", func=self.get_cell_params, dependencies=[dget(self.cell, "h")])

      dset(self, "stress", depend_value(name="stress", func=self.get_stress, dependencies=[dget(self.beads, "kstress"), dget(self.forces, "vir"), dget(self.cell, "V")]))
      self.property_dict["stress_md.xx"] = depend_value(name="scl_xx", dependencies=[dget(self, "stress")], func=(lambda : self.stress[0,0]) ) 
      
      self.add_property(prop_name="pressure_md", dep_name="press", func=self.get_press, dependencies=[dget(self, "stress")])
      self.add_property(prop_name="kinetic_cv", dep_name="kin_cv", func=self.get_kincv, dependencies=[dget(self.beads, "q"), dget(self.forces, "f"), dget(self.ensemble, "temp")])

      dset(self, "kstress_cv", depend_value(name="kstress_cv", func=self.get_kstresscv, dependencies=[dget(self.beads,"q"),dget(self.forces,"f"),dget(self.ensemble,"temp")]))
      dset(self, "stress_cv", depend_value(name="stress_cv", func=self.get_stresscv, dependencies=[dget(self,"kstress_cv"),dget(self.forces,"vir"), dget(self.cell, "V")]))
      self.property_dict["stress_cv.xx"] = depend_value(name="scv_xx", dependencies=[dget(self, "stress_cv")], func=(lambda : self.stress_cv[0,0]) ) 

      self.add_property(prop_name="pressure_cv", dep_name="press_cv", func=self.get_presscv, dependencies=[dget(self, "stress_cv")])
      
      self.dbeads = simul.beads.copy()
      self.dforces = ForceBeads()
      self.dforces.bind(self.dbeads, self.simul.cell,  self.simul._forcemodel)
      self.add_property(prop_name="kinetic_yamamoto", dep_name="kin_yama", func=self.get_kinyama, dependencies=[dget(self.beads, "q"), dget(self.ensemble, "temp")])
      
   def add_property(self, prop_name, dep_name, func, dependencies=None):
      """Adds a property to the property list.

      Args:
         prop_name: A string giving the keyword that will identify the property 
            in the output list.
         dep_name: A string giving the name of the attribute created in the 
            Properties object.
         func: The function used to compute the property.
         dependencies: A list of depend objects on which the property will 
            depend upon.
      """

      dset(self, dep_name, depend_value(name=dep_name, func=func, dependencies=dependencies))
      self.property_dict[prop_name] = dget(self, dep_name)

   def get_kin(self):
      """Calculates the classical kinetic energy estimator."""

      return self.beads.kin/self.beads.nbeads

   def get_time(self):
      """Calculates the elapsed simulation time."""

      return (1 + self.simul.step)*self.ensemble.dt

   def get_step(self):
      """Return the simulation step."""

      return (1 + self.simul.step)

   def __getitem__(self,key):
      """Retrieves the item given by key.

      Args:
         key: A string contained in property_dict.

      Returns:
         The property labelled by the keyword key.
      """

      return self.property_dict[key].get()

   def get_pot(self):
      """Calculates the potential energy estimator."""

      return self.forces.pot/self.beads.nbeads

   def get_temp(self):
      """Calculates the classical kinetic temperature estimator.

      Note that in the case that the centre of mass constraint there will be
      one less degree of freedom than without, so this has to be taken into
      account when calculating the kinetic temperature.
      """
      if self.ensemble.fixcom: mdof=3 
      else: mdof=0
      return self.beads.kin/(0.5*Constants.kb*(3*self.beads.natoms*self.beads.nbeads - mdof)*self.beads.nbeads)

   def get_econs(self):
      """Calculates the conserved quantity estimator."""

      return self.ensemble.econs/(self.beads.nbeads*self.beads.natoms)

   def get_stress(self):
      """Calculates the classical kinetic energy estimator."""

      return (self.forces.vir + self.beads.kstress)/self.cell.V

   def get_press(self):
      """Calculates the classical pressure estimator."""

      return np.trace(self.stress)/3.0

   def get_stresscv(self):
      """Calculates the quantum central virial stress tensor estimator."""

      return (self.forces.vir + self.kstress_cv)/self.cell.V                  

   def get_presscv(self):
      """Calculates the quantum central virial pressure estimator."""

      return np.trace(self.stress_cv)/3.0
   
   def get_kincv(self):        
      """Calculates the quantum central virial kinetic energy estimator."""

      kcv=0.0
      for b in range(self.beads.nbeads):
         kcv += np.dot(depstrip(self.beads.q[b]) - depstrip(self.beads.qc), depstrip(self.forces.f[b]))
      kcv *= -0.5/self.beads.nbeads
      kcv += 0.5*Constants.kb*self.ensemble.temp*(3*self.beads.natoms) 
      return kcv

   def get_kstresscv(self):        
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

   def get_kinyama(self):              
      """Calculates the quantum scaled coordinate kinetic energy estimator.

      Uses a finite difference method to calculate the kinetic energy estimator
      without requiring the forces as for the centroid virial estimator.
      """
      
      dbeta = abs(self.fd_delta)
      
      v0 = self.pot
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
         
   def get_cell_params(self):
      """Returns a list of the cell box lengths and the angles between them.

      Returns:
         A list of the form [a, b, c, alpha, beta, gamma].
      """

      a, b, c, alpha, beta, gamma = h2abc(self.cell.h)
      return [a, b, c, alpha, beta, gamma]


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
