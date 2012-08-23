"""Contains the classes that deal with the different dynamics required in
different types of ensembles.

Holds the algorithms required for normal mode propagators, and the objects to
do the constant temperature and pressure algorithms. Also calculates the
appropriate conserved energy quantity for the ensemble of choice.

Classes:
   Ensemble: Base ensemble class with generic methods and attributes.
   NVEEnsemble: Deals with constant energy dynamics.
   NVTEnsemble: Deals with constant temperature dynamics.
   NPTEnsemble: Deals with constant pressure dynamics.
"""

__all__ = ['Ensemble', 'NVEEnsemble', 'NVTEnsemble', 'NPTEnsemble']

import numpy as np
import math
from utils.depend import *
from utils import units
from thermostats import *
from barostats import *
from inputs.thermostats import InputThermo
from inputs.barostats import InputBaro
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
         nm_freqs:
      """
#TODO add nm_freqs to the doc string

      dset(self, "econs", depend_value(name='econs', func=self.get_econs) )
      dset(self, "temp",  depend_value(name='temp',  value=temp))
      dset(self, "dt",    depend_value(name='dt',    value=dt))


   def bind(self, beads, nm, cell, bforce, prng):
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
      self.nm = nm

      # n times the temperature
      dset(self,"ntemp", depend_value(name='ntemp',func=self.get_ntemp,dependencies=[dget(self,"temp")]))

      # dependencies of the conserved quantity
      dget(self,"econs").add_dependency(dget(self.beads, "kin"))
      dget(self,"econs").add_dependency(dget(self.forces, "pot"))
      dget(self,"econs").add_dependency(dget(self.beads, "vpath"))


   def get_ntemp(self):
      """Returns the PI simulation temperature (P times the physical T)."""

      return self.temp*self.beads.nbeads


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

      return self.nm.kin + self.beads.vpath*self.nm.omegan2 + self.forces.pot


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

      self.nm.qnm[0,:] += depstrip(self.nm.pnm)[0,:]/depstrip(self.beads.m3)[0]*self.dt

   def step(self):
      """Does one simulation time step."""

      self.ptime = -time.time()
      self.pstep()
      self.ptime += time.time()

      self.qtime = -time.time()
      self.qcstep()
      self.nm.free_qstep()
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

   def bind(self, beads, nm, cell, bforce, prng):
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

      super(NVTEnsemble,self).bind(beads, nm, cell, bforce, prng)
      ndof = None
      if self.fixcom:
         ndof = 3*(self.beads.natoms-1)

      if isinstance(self.thermostat,ThermoNMGLE) or isinstance(self.thermostat,ThermoNMGLEG) or isinstance(self.thermostat,ThermoPILE_L) or isinstance(self.thermostat,ThermoPILE_G):
         self.thermostat.bind(nm=self.nm,prng=prng,ndof=ndof )
      else:
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
      self.nm.free_qstep()
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


class NPTEnsemble(NVTEnsemble):
   """Ensemble object for constant pressure simulations.

   Has the relevant conserved quantity and normal mode propagator for the
   constant pressure ensemble. Contains a thermostat object containing the
   algorithms to keep the temperature constant, and a barostat to keep the
   pressure constant.

   Attributes:
      barostat: A barostat object to keep the pressure constant.

   Depend objects:
      econs: Conserved energy quantity. Depends on the bead and cell kinetic
         and potential energy, the spring potential energy, the heat
         transferred to the beads and cell thermostat, the temperature and
         the cell volume.
      pext: External pressure.
      sext: External stress tensor.
   """

   def __init__(self, dt, temp, pext=None, sext=None, thermostat=None, barostat=None, nm_freqs=None, fixcom=False):
      """Initialises NPTEnsemble.

      Args:
         dt: The simulation timestep.
         temp: The system temperature.
         pext: An optional float giving the external pressure.
         sext: An optional array giving the external stress tensor.
         thermostat: A thermostat object to keep the temperature constant.
            Defaults to Thermostat().
         barostat: A barostat object to keep the temperature constant.
            Defaults to Barostat().
         fixcom: An optional boolean which decides whether the centre of mass
            motion will be constrained or not. Defaults to False.
      """

      super(NPTEnsemble,self).__init__(dt, temp, thermostat, nm_freqs=nm_freqs, fixcom=fixcom)
      if barostat == None:
         self.barostat = Barostat()
      else:
         self.barostat = barostat

      sync_ext=synchronizer()
      dset(self,"sext",
         depend_array(name='sext', value=np.zeros((3,3)), synchro=sync_ext,
            func={"pext" : self.p2s}))
      dset(self,"pext",
         depend_value(name='pext', value=0.0, synchro=sync_ext,
            func={"sext" : self.s2p}))
      if pext is not None:
         dset(self,"pext",depend_value(name="pext", value=pext))
         deppipe(self, "pext", self.barostat, "pext")
         deppipe(self, "sext", self.barostat, "sext")
      elif sext is not None:
         dset(self,"sext",depend_value(name="sext", value=pext))
         deppipe(self, "sext", self.barostat, "sext")
         deppipe(self, "pext", self.barostat, "pext")
      else:
         raise TypeError("You must provide either the pressure or stress")

   def s2p(self):
      """Converts the external stress to the external pressure."""

      return np.trace(self.sext)/3.0

   def p2s(self):
      """Converts the external pressure to an isotropic external stress."""

      return self.pext*np.identity(3)

   def bind(self, beads, nm, cell, bforce, prng):
      """Binds beads, cell, bforce and prng to the ensemble.

      This takes a beads object, a cell object, a forcefield object and a
      random number generator object and makes them members of the ensemble.
      It also then creates the objects that will hold the data needed in the
      ensemble algorithms and the dependency network. Also note that the cell
      thermostat timesteps and temperatures are defined relative to the system
      temperature, and the the thermostat temperatures are held at the
      higher simulation temperature, as is appropriate.

      Args:
         beads: The beads object from whcih the bead positions are taken.
         cell: The cell object from which the system box is taken.
         bforce: The forcefield object from which the force and virial are
            taken.
         prng: The random number generator object which controls random number
            generation.
      """

      super(NPTEnsemble,self).bind(beads, cell, nm, bforce, prng)
      self.barostat.bind(beads, cell, bforce)

      deppipe(self,"ntemp", self.barostat,"temp")
      deppipe(self,"ntemp", self.barostat.thermostat,"temp")
      deppipe(self, "dt", self.barostat, "dt")

      dget(self,"econs").add_dependency(dget(self.barostat.thermostat, "ethermo"))
      dget(self,"econs").add_dependency(dget(self.barostat, "pot"))
      dget(self,"econs").add_dependency(dget(self.thermostat, "temp"))
      dget(self,"econs").add_dependency(dget(self.cell, "kin"))
      dget(self,"econs").add_dependency(dget(self.cell, "V"))

   def get_econs(self):
      """Calculates the conserved energy quantity for the constant pressure
      ensemble.
      """

      return NVTEnsemble.get_econs(self) + self.barostat.thermostat.ethermo + self.barostat.pot + self.cell.kin - 2.0*units.Constants.kb*self.thermostat.temp*math.log(self.cell.V)

   def step(self):
      """NPT time step.

      Note that the barostat only propagates the centroid coordinates. If this
      approximation is made a centroid virial pressure and stress estimator can
      be defined, so this gives the best statistical convergence. This is
      allowed as the normal mode propagation is approximately unaffected
      by volume fluctuations as long as the system box is much larger than
      the radius of gyration of the ring polymers.
      """

      self.ttime = -time.time()
      self.thermostat.step()
      self.barostat.thermostat.step()
      self.rmcom()
      self.ttime += time.time()

      self.ptime = -time.time()
      self.barostat.pstep()
      self.ptime += time.time()

      self.qtime = -time.time()
      self.barostat.qcstep()
      self.qstep()
      self.qtime += time.time()

      self.ptime -= time.time()
      self.barostat.pstep()
      self.ptime += time.time()

      self.ttime -= time.time()
      self.barostat.thermostat.step()
      self.thermostat.step()
      self.rmcom()
      self.ttime += time.time()
