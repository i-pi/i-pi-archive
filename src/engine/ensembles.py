"""Contains the classes that deal with the different dynamics required in 
different types of ensembles.

Holds the algorithms required for normal mode propagators, and the objects to
do the constant temperature algorithm. Also calculates the 
appropriate conserved energy quantity for the ensemble of choice.

Classes:
   Ensemble: Base ensemble class with generic methods and attributes.
   NVEEnsemble: Deals with constant energy dynamics.
   NVTEnsemble: Deals with constant temperature dynamics.
"""

__all__ = ['Ensemble', 'NVEEnsemble', 'NVTEnsemble' ]

import numpy as np
import math 
from utils.depend import *
from utils import units
from thermostats import *
from inputs.thermostats import InputThermo
import time

class Ensemble(dobject): 
   """Base (do-nothing) ensemble class.

      Gives the standard methods and attributes needed in all the 
      ensemble classes.

      Attributes:
         beads: A beads object giving the atoms positions.
         cell: A cell object giving the system box.
         forces: A forces object giving the virial and the forces acting on
            each bead.
         prng: A random number generator object.

      Depend objects:
         econs: The conserved energy quantity appropriate to the given 
            ensemble. Depends on the various energy terms which make it up,
            which are different depending on the ensemble.
         temp: The system temperature.
         dt: The timestep for the algorithms.
         ntemp: The simulation temperature. Will be nbeads times higher than
            the system temperature as PIMD calculations are done at this 
            effective classical temperature.
         omegan: The effective vibrational frequency for the interaction 
            between the replicas. Depends on the simulation temperature.
         omegan2: omegan**2.
         omegak: The normal mode frequencies for the free ring polymer.
            Depends on omegan.
         prop_pq: An array holding the exact normal mode propagator for the
            free ring polymer, using mass scaled coordinates. 
            See J. Chem. Phys. 133, 124101 (2010). Depends on the bead masses
            and the timestep.
      """

   def __init__(self, dt, temp):
      """Initialises Ensemble.

      Args:
         dt: The timestep of the simulation algorithms.
         temp: The temperature.
      """

      dset(self, "econs", depend_value(name='econs', func=self.get_econs) )
      dset(self, "temp",  depend_value(name='temp',  value=temp))       
      dset(self, "dt",    depend_value(name='dt',    value=dt))
      
   def bind(self, beads, cell, bforce, prng):
      """Binds beads, cell, bforce and prng to the ensemble.

      This takes a beads object, a cell object, a forcefield object and a 
      random number generator object and makes them members of the ensemble.
      It also then creates the objects that will hold the data needed in the
      ensemble algorithms and the dependency network. Note that the conserved
      quantity is defined in the init, but as each ensemble has a different
      conserved quantity the dependencies are defined in bind.

      Args:
         beads: The beads object from whcih the bead positions are taken.
         cell: The cell object from which the system box is taken.
         bforce: The forcefield object from which the force and virial are
            taken.
         prng: The random number generator object which controls random number
            generation.
      """

      # store local references to the different bits of the simulation
      self.beads = beads
      self.cell = cell
      self.forces = bforce
      self.prng = prng
       
      # n times the temperature 
      dset(self,"ntemp", depend_value(name='ntemp',func=self.get_ntemp,dependencies=[dget(self,"temp")]))
      
      # dependencies of the conserved quantity
      dget(self,"econs").add_dependency(dget(self.beads, "kin"))
      dget(self,"econs").add_dependency(dget(self.forces, "pot"))
      dget(self,"econs").add_dependency(dget(self.beads, "vpath"))
      
      # create path related properties
      dset(self,"omegan",
         depend_value(name='omegan',func=self.get_omegan, 
            dependencies=[dget(self,"ntemp")]) )
      dset(self,"omegan2",
         depend_value(name='omegan2',func=self.get_omegan2, 
            dependencies=[dget(self,"omegan")]) )
      dset(self,"omegak",
         depend_array(name='omegak',value=np.zeros(self.beads.nbeads,float),
            func=self.get_omegak, dependencies=[dget(self,"omegan")]) )
      dset(self,"prop_pq",
         depend_array(name='prop_pq',value=np.zeros((self.beads.nbeads,2,2)),
            func=self.get_prop_pq, 
               dependencies=[dget(self,"omegak"), dget(self,"dt")]) )

   def get_ntemp(self):
      """Returns the PI simulation temperature (P times the physical T)."""

      return self.temp*self.beads.nbeads

   def get_omegan(self):
      """Returns the effective vibrational frequency for the interaction 
      between replicas.
      """

      return self.ntemp*units.Constants.kb/units.Constants.hbar

   def get_omegan2(self):
      """Returns omegan**2."""

      return self.omegan**2

   def get_omegak(self):
      """Gets the normal mode frequencies.

      Returns:
         A list of the normal mode frequencies for the free ring polymer.
         The first element is the centroid frequency (0.0).
      """

      return 2*self.omegan*np.array([math.sin(k*math.pi/self.beads.nbeads) for k in range(self.beads.nbeads)])

   def get_prop_pq(self): 
      """Gets the normal mode propagator matrix.

      Note the special treatment for the centroid normal mode, which is
      propagated using the standard velocity Verlet algorithm as required.
      Note that both the normal mode positions and momenta are propagated
      using this matrix.

      Returns:
         An array of the form (nbeads, 2, 2). Each 2*2 array prop_pq[i,:,:] 
         gives the exact propagator for the i-th normal mode of the 
         ring polymer.
      """

      pqk = np.zeros((self.beads.nbeads,2,2), float)
      pqk[0] = np.array([[1,0], [self.dt,1]])
      for b in range(1, self.beads.nbeads):
         dtomegak = self.omegak[b]*self.dt
         c = math.cos(dtomegak)
         s = math.sin(dtomegak)
         pqk[b,0,0] = c
         pqk[b,1,1] = c
         pqk[b,0,1] = -s*self.omegak[b]
         pqk[b,1,0] = s/self.omegak[b]
      return pqk

         
   def pstep(self): 
      """Dummy momenta propagator which does nothing."""
      
      pass

   def qcstep(self): 
      """Dummy centroid position propagator which does nothing."""

      pass

   def step(self): 
      """Dummy simulation time step which does nothing."""

      pass

   def get_econs(self):
      """Calculates the conserved energy quantity for constant energy 
      ensembles.
      """

      return self.beads.kin + self.beads.vpath*self.omegan2 + self.forces.pot
      
      
class NVEEnsemble(Ensemble):
   """Ensemble object for constant energy simulations.

   Has the relevant conserved quantity and normal mode propagator for the 
   constant energy ensemble. Note that a temperature of some kind must be 
   defined so that the spring potential can be calculated.

   Attributes:
      fixcom: A boolean which decides whether the centre of mass 
         motion will be constrained or not.
   
   Depend objects:
      econs: Conserved energy quantity. Depends on the bead kinetic and 
         potential energy, and the spring potential energy.
   """

   def __init__(self, dt, temp, fixcom=False):
      """Initialises NVEEnsemble.

      Args:
         dt: The simulation timestep.
         temp: The system temperature.
         fixcom: An optional boolean which decides whether the centre of mass
            motion will be constrained or not. Defaults to False.
      """

      super(NVEEnsemble,self).__init__(dt=dt,temp=temp)
      self.fixcom = fixcom

   def rmcom(self):
      """This removes the centre of mass contribution to the kinetic energy.

      Calculates the centre of mass momenta, then removes the mass weighted
      contribution from each atom. If the ensemble defines a thermostat, then
      the contribution to the conserved quantity due to this subtraction is 
      added to the thermostat heat energy, as it is assumed that the centre of
      mass motion is due to the thermostat. 

      If there is a choice of thermostats, the thermostat 
      connected to the centroid is chosen.
      """

      if (self.fixcom):
         pcom = np.zeros(3,float);
         
         p = depstrip(self.beads.p)
         na3 = self.beads.natoms*3
         nb = self.beads.nbeads
         m = depstrip(self.beads.m3)[:,0:na3:3]
         for i in range(3):
            pcom[i] = p[:,i:na3:3].sum()

         if hasattr(self,"thermostat"):
            if hasattr(self.thermostat, "_thermos"):
               self.thermostat._thermos[0].ethermo += np.dot(pcom,pcom)/(2.0*self.beads[0].M*nb)
            else:
               self.thermostat.ethermo += np.dot(pcom,pcom)/(2.0*self.beads[0].M*nb)

         # subtracts COM _velocity_
         pcom *= 1.0/(nb*self.beads[0].M)
         for i in range(3):
            self.beads.p[:,i:na3:3] -= m*pcom[i]

   def pstep(self): 
      """Velocity Verlet momenta propagator."""

      self.beads.p += depstrip(self.forces.f)*(self.dt*0.5)
   
   def qcstep(self):
      """Velocity Verlet centroid position propagator."""
      self.beads.qnm[0,:] += depstrip(self.beads.pnm)[0,:]/depstrip(self.beads.m3)[0]*self.dt

   def qstep(self):
      """Exact normal mode propagator for the free ring polymer.

      Note that the propagator works in mass scaled coordinates, so that the
      propagator matrix can be determined independently from the particular
      atom masses, and so the same propagator will work for all the atoms in 
      the system. All the ring polymers are propagated at the same time by a
      matrix multiplication.

      Also note that the centroid coordinate is propagated in qcstep, so is
      not altered here.
      """

      if self.beads.nbeads == 1:
         pass
      else:
         pq = np.zeros((2,self.beads.natoms*3),float)
         sm = depstrip(self.beads.sm3[0])
         for k in range(1,self.beads.nbeads):
            pq[0,:] = depstrip(self.beads.pnm[k])/sm
            pq[1,:] = depstrip(self.beads.qnm[k])*sm
            pq = np.dot(self.prop_pq[k],pq)
            self.beads.qnm[k] = pq[1,:]/sm         
            self.beads.pnm[k] = pq[0,:]*sm                    
      
   def step(self): 
      """Does one simulation time step."""

      self.ptime = -time.time()
      self.pstep()
      self.ptime += time.time()

      self.qtime = -time.time()
      self.qcstep()
      self.qstep()
      self.qtime += time.time()

      self.ptime -= time.time()
      self.pstep()
      self.ptime += time.time()

      self.ttime = -time.time()
      self.rmcom()
      self.ttime += time.time()
      

class NVTEnsemble(NVEEnsemble):
   """Ensemble object for constant temperature simulations.

   Has the relevant conserved quantity and normal mode propagator for the 
   constant temperature ensemble. Contains a thermostat object containing the
   algorithms to keep the temperature constant.

   Attributes:
      thermostat: A thermostat object to keep the temperature constant.
   
   Depend objects:
      econs: Conserved energy quantity. Depends on the bead kinetic and 
         potential energy, the spring potential energy and the heat 
         transferred to the thermostat.
   """

   def __init__(self, dt, temp, thermostat=None, fixcom=False):
      """Initialises NVTEnsemble.

      Args:
         dt: The simulation timestep.
         temp: The system temperature.
         thermostat: A thermostat object to keep the temperature constant.
            Defaults to Thermostat()
         fixcom: An optional boolean which decides whether the centre of mass
            motion will be constrained or not. Defaults to False.
      """

      super(NVTEnsemble,self).__init__(dt=dt,temp=temp, fixcom=fixcom)

      if thermostat is None:
         self.thermostat = Thermostat()
      else:
         self.thermostat = thermostat

   def bind(self, beads, cell, bforce, prng):
      """Binds beads, cell, bforce and prng to the ensemble.

      This takes a beads object, a cell object, a forcefield object and a 
      random number generator object and makes them members of the ensemble.
      It also then creates the objects that will hold the data needed in the
      ensemble algorithms and the dependency network. Also note that the 
      thermostat timestep and temperature are defined relative to the system
      temperature, and the the thermostat temperature is held at the 
      higher simulation temperature, as is appropriate.

      Args:
         beads: The beads object from whcih the bead positions are taken.
         cell: The cell object from which the system box is taken.
         bforce: The forcefield object from which the force and virial are
            taken.
         prng: The random number generator object which controls random number
            generation.
      """

      super(NVTEnsemble,self).bind(beads, cell, bforce, prng)
      ndof = None 
      if self.fixcom:
         ndof = 3*(self.beads.natoms-1)
         
      self.thermostat.bind(beads=self.beads,prng=prng,ndof=ndof )

      deppipe(self,"ntemp", self.thermostat,"temp")
      deppipe(self,"dt", self.thermostat, "dt")

      dget(self,"econs").add_dependency(dget(self.thermostat, "ethermo"))      
      
   def step(self): 
      """Does one simulation time step."""

      self.ttime = -time.time()
      self.thermostat.step()
      self.rmcom()
      self.ttime += time.time()

      self.ptime = -time.time()
      self.pstep()
      self.ptime += time.time()
      
      self.qtime = -time.time()
      self.qcstep()
      self.qstep()
      self.qtime += time.time()

      self.ptime -= time.time()
      self.pstep()
      self.ptime += time.time()

      self.ttime -= time.time()
      self.thermostat.step()
      self.rmcom()
      self.ttime += time.time()

   def get_econs(self):
      """Calculates the conserved energy quantity for constant temperature 
      ensemble.
      """

      return NVEEnsemble.get_econs(self) + self.thermostat.ethermo 
