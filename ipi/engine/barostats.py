"""Classes that deal with constant pressure simulations.

Contains the algorithms which propagate the position and momenta steps in the
constant pressure ensemble. Holds the properties directly related to
these ensembles, such as the internal and external pressure and stress.

Note that this file also contains a 'BaroMHT' class, that follows more closely the
Martyna, Hughes, Tuckerman implementation of a PIMD barostat. However it is so
close to the BZP implementation that we disabled it for the sake of simplicity.
The original reference is:
G. Martyna, A. Hughes and M. Tuckerman, J. Chem. Phys., 110, 3275.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

from ipi.utils.depend import *
from ipi.utils.units import *
from ipi.utils.mathtools import eigensystem_ut3x3, invert_ut3x3, exp_ut3x3, det_ut3x3, matrix_exp
from ipi.inputs.thermostats import InputThermo
from ipi.engine.thermostats import Thermostat
from ipi.engine.cell import Cell


__all__ = ['Barostat', 'BaroBZP', 'BaroRGB']


class Barostat(dobject):
   """Base barostat class.

   Gives the standard methods and attributes needed in all the barostat classes.

   Attributes:
      beads: A beads object giving the atoms positions
      cell: A cell object giving the system box.
      forces: A forces object giving the virial and the forces acting on
         each bead.
      nm: An object to do the normal mode transformation.
      thermostat: A thermostat coupled to the barostat degrees of freedom.
      mdof: The number of atomic degrees of freedom

   Depend objects:
      dt: The time step used in the algorithms. Depends on the simulation dt.
      temp: The (classical) simulation temperature. Higher than the physical
         temperature by a factor of the number of beads.
      tau: The timescale associated with the piston
      pext: The external pressure
      stressext: The external stress
      ebaro: The conserved quantity associated with the barostat.
      pot: The potential energy associated with the barostat.
      kstress: The system kinetic stress tensor.
      stress: The system stress tensor.
      press: The system pressure.
   """

   def __init__(self, dt=None, temp=None, tau=None, ebaro=None, thermostat=None):
      """Initialises base barostat class.

      Note that the external stress and the external pressure are synchronized.
      This makes most sense going from the stress to the pressure, but if you
      must go in the other direction the stress is assumed to be isotropic.

      Args:
         dt: Optional float giving the time step for the algorithms. Defaults
            to the simulation dt.
         temp: Optional float giving the temperature for the thermostat.
            Defaults to the simulation temp.
         tau: Optional float giving the time scale associated with the barostat.
         ebaro: Optional float giving the conserved quantity already stored
            in the barostat initially. Used on restart.
         thermostat: The thermostat connected to the barostat degree of freedom.
      """

      dset(self,"dt",depend_value(name='dt'))
      if not dt is None:
         self.dt = dt
      else: self.dt = 1.0

      dset(self, "temp", depend_value(name="temp"))
      if not temp is None:
         self.temp = temp
      else: self.temp = 1.0

      dset(self,"tau",depend_value(name='tau'))
      if not tau is None:
         self.tau = tau
      else: self.tau = 1.0

      dset(self,"ebaro",depend_value(name='ebaro'))
      if not ebaro is None:
         self.ebaro = ebaro
      else: self.ebaro = 0.0

      if thermostat is None:
         thermostat = Thermostat()
      self.thermostat = thermostat

      # pipes timestep and temperature to the thermostat
      deppipe(self, "dt", self.thermostat, "dt")
      deppipe(self, "temp", self.thermostat, "temp")
      dset(self, "pext", depend_value(name='pext', value=-1.0))
      dset(self, "stressext", depend_array(name='stressext', value=-np.ones((3,3), float)))


   def bind(self, beads, nm, cell, forces, bias=None, prng=None, fixdof=None):
      """Binds beads, cell and forces to the barostat.

      This takes a beads object, a cell object and a forcefield object and
      makes them members of the barostat. It also then creates the objects that
      will hold the data needed in the barostat algorithms and the dependency
      network.

      Args:
         beads: The beads object from which the bead positions are taken.
         nm: The normal modes propagator object
         cell: The cell object from which the system box is taken.
         forces: The forcefield object from which the force and virial are
            taken.
         prng: The parent PRNG to bind the thermostat to
         fixdof: The number of blocked degrees of freedom.
      """

      self.beads = beads
      self.cell = cell
      self.forces = forces
      self.bias = bias
      self.nm = nm

      dset(self,"kstress",
         depend_value(name='kstress', func=self.get_kstress,
            dependencies=[ dget(beads,"q"), dget(beads,"qc"), dget(beads,"pc"), dget(forces,"f") ]))
      dset(self,"stress",
         depend_value(name='stress', func=self.get_stress,
            dependencies=[ dget(self,"kstress"), dget(cell,"V"), dget(forces,"vir") ]))
      if bias != None:
         dget(self,"kstress").add_dependency(dget(bias,"f"))
         dget(self,"stress").add_dependency(dget(bias,"vir"))

      if fixdof is None:
         self.mdof = float(self.beads.natoms)*3.0
      else:
         self.mdof = float(self.beads.natoms)*3.0 - float(fixdof)


   def get_kstress(self):
      """Calculates the quantum centroid virial kinetic stress tensor
      estimator.
      """

      kst = np.zeros((3,3),float)
      q = depstrip(self.beads.q)
      qc = depstrip(self.beads.qc)
      pc = depstrip(self.beads.pc)
      m = depstrip(self.beads.m)
      na3 = 3*self.beads.natoms
      fall = depstrip(self.forces.f)
      if self.bias == None: ball=fall*0.00
      else: ball = depstrip(self.bias.f)

      for b in range(self.beads.nbeads):
         for i in range(3):
            for j in range(i,3):
               kst[i,j] -= np.dot(q[b,i:na3:3] - qc[i:na3:3],
                  fall[b,j:na3:3]+ball[b,j:na3:3])

      # NOTE: In order to have a well-defined conserved quantity, the Nf kT term in the
      # diagonal stress estimator must be taken from the centroid kinetic energy.
      for i in range(3):
         kst[i,i] += np.dot(pc[i:na3:3],pc[i:na3:3]/m) *self.beads.nbeads

      return kst

   def get_stress(self):
      """Calculates the internal stress tensor."""

      bvir = np.zeros((3,3),float)
      if self.bias != None: bvir[:]=self.bias.vir
      return (self.kstress + self.forces.vir + bvir)/self.cell.V

   def pstep(self):
      """Dummy momenta propagator step."""

      pass

   def qcstep(self):
      """Dummy centroid position propagator step."""

      pass


class BaroBZP(Barostat):
   """Bussi-Zykova-Parrinello barostat class.

   Just extends the standard class adding finite-dt propagators for the
   barostat velocities, positions, piston.

   Generates dynamics with a stochastic barostat. Implementation details:
   Ceriotti, More, Manolopoulos, Comp. Phys. Comm. 185, 1019, (2013)

   Depend objects:
      p: The momentum associated with the volume degree of freedom.
      m: The mass associated with the volume degree of freedom.
   """

   def __init__(self, dt=None, temp=None, tau=None, ebaro=None, thermostat=None, pext=None, p=None):
      """Initializes BZP barostat.

      Args:
         dt: Optional float giving the time step for the algorithms. Defaults
            to the simulation dt.
         temp: Optional float giving the temperature for the thermostat.
            Defaults to the simulation temp.
         pext: Optional float giving the external pressure.
         tau: Optional float giving the time scale associated with the barostat.
         ebaro: Optional float giving the conserved quantity already stored
            in the barostat initially. Used on restart.
         thermostat: The thermostat connected to the barostat degree of freedom.
         p: Optional initial volume conjugate momentum. Defaults to 0.
      """


      super(BaroBZP, self).__init__(dt, temp, tau, ebaro, thermostat)

      dset(self,"p", depend_array(name='p', value=np.atleast_1d(0.0)))

      if not p is None:
         self.p = np.asarray([p])
      else:
         self.p = 0.0
      
      if not pext is None:
         self.pext = pext
      else:
         self.pext = -1.0

   def bind(self, beads, nm, cell, forces, bias=None, prng=None, fixdof=None):
      """Binds beads, cell and forces to the barostat.

      This takes a beads object, a cell object and a forcefield object and
      makes them members of the barostat. It also then creates the objects that
      will hold the data needed in the barostat algorithms and the dependency
      network.

      Args:
         beads: The beads object from which the bead positions are taken.
         nm: The normal modes propagator object
         cell: The cell object from which the system box is taken.
         forces: The forcefield object from which the force and virial are
            taken.
         prng: The parent PRNG to bind the thermostat to
         fixdof: The number of blocked degrees of freedom.
      """

      super(BaroBZP, self).bind(beads, nm, cell, forces, bias, prng, fixdof)

      # obtain the thermostat mass from the given time constant
      # note that the barostat temperature is nbeads times the physical T
      dset(self,"m", depend_array(name='m', value=np.atleast_1d(0.0),
         func=(lambda:np.asarray([self.tau**2*3*self.beads.natoms*Constants.kb*self.temp])),
            dependencies=[ dget(self,"tau"), dget(self,"temp") ] ))

      # binds the thermostat to the piston degrees of freedom      
      self.thermostat.bind(pm=[ self.p, self.m ], prng=prng)

      # barostat elastic energy
      dset(self,"pot",
         depend_value(name='pot', func=self.get_pot,
            dependencies=[ dget(cell,"V"), dget(self,"pext") ]))

      dset(self,"kin",depend_value(name='kin',
         func=(lambda:0.5*self.p[0]**2/self.m[0]),
            dependencies= [dget(self,"p"), dget(self,"m")] ) )

      # the barostat energy must be computed from bits & pieces (overwrite the default)
      dset(self, "ebaro", depend_value(name='ebaro', func=self.get_ebaro,
         dependencies=[ dget(self, "kin"), dget(self, "pot"),
            dget(self.cell, "V"), dget(self, "temp"),
               dget(self.thermostat,"ethermo")] ))

   def get_pot(self):
      """Calculates the elastic strain energy of the cell."""

      # NOTE: since there are nbeads replicas of the unit cell, the enthalpy contains a nbeads factor
      return self.cell.V*self.pext*self.beads.nbeads

   def get_ebaro(self):
      """Calculates the barostat conserved quantity."""

      return self.thermostat.ethermo + self.kin + self.pot - np.log(self.cell.V)*Constants.kb*self.temp

   def pstep(self):
      """Propagates the momenta for half a time step."""

      dthalf = self.dt*0.5
      dthalf2 = dthalf**2
      dthalf3 = dthalf**3/3.0

      press = np.trace(self.stress)/3.0
      # This differs from the BZP thermostat in that it uses just one kT in the propagator.
      # This leads to an ensemble equaivalent to Martyna-Hughes-Tuckermann for both fixed and moving COM
      # Anyway, it is a small correction so whatever.
      self.p += dthalf*3.0*( self.cell.V* ( press - self.beads.nbeads*self.pext ) +
                Constants.kb*self.temp )

      fc = np.sum(depstrip(self.forces.f),0)/self.beads.nbeads
      if self.bias != None: fc += np.sum(depstrip(self.bias.f),0)/self.beads.nbeads
      m = depstrip(self.beads.m3)[0]
      pc = depstrip(self.beads.pc)

      # I am not 100% sure, but these higher-order terms come from integrating the pressure virial term,
      # so they should need to be multiplied by nbeads to be consistent with the equations of motion in the PI context
      # again, these are tiny tiny terms so whatever.
      self.p += (dthalf2*np.dot(pc,fc/m) + dthalf3*np.dot(fc,fc/m)) * self.beads.nbeads

      self.beads.p += depstrip(self.forces.f)*dthalf
      if self.bias != None: self.beads.p +=depstrip(self.bias.f)*dthalf

   def qcstep(self):
      """Propagates the centroid position and momentum and the volume."""

      v = self.p[0]/self.m[0]
      expq, expp = (np.exp(v*self.dt), np.exp(-v*self.dt))

      m = depstrip(self.beads.m3)[0]

      self.nm.qnm[0,:] *= expq
      self.nm.qnm[0,:] += ((expq-expp)/(2.0*v))* (depstrip(self.nm.pnm)[0,:]/m)
      self.nm.pnm[0,:] *= expp

      self.cell.h *= expq


class BaroRGB(Barostat):
   """Raiteri-Gale-Bussi constant stress barostat class (JPCM 23, 334213, 2011).

      Just extends the standard class adding finite-dt propagators for the barostat
      velocities, positions, piston.

      Depend objects:
      p: The momentum matrix associated with the cell degrees of freedom.
      m: The mass associated with the cell degree of freedom.
      """

   def __init__(self, dt=None, temp=None, tau=None, ebaro=None, thermostat=None, stressext=None, h0=None, p=None):
      """Initializes BZP barostat.

         Args:
         dt: Optional float giving the time step for the algorithms. Defaults
         to the simulation dt.
         temp: Optional float giving the temperature for the thermostat.
         Defaults to the simulation temp.
         stressext: Optional float giving the external pressure.
         tau: Optional float giving the time scale associated with the barostat.
         ebaro: Optional float giving the conserved quantity already stored
         in the barostat initially. Used on restart.
         thermostat: The thermostat connected to the barostat degree of freedom.
         p: Optional initial volume conjugate momentum. Defaults to 0.
         """


      super(BaroRGB, self).__init__(dt, temp, tau, ebaro, thermostat)

      # non-zero elements of the cell momentum are only
      # pxx pyy pzz pxy pxz pyz, but we want to access it either as a
      # 6-vector or as a 3x3 upper triangular tensor.
      # we use a synchronizer to achieve that

      sync_baro = synchronizer()
      dset(self,"p6", depend_array(name='p6', value=np.zeros(6,float),
          synchro=sync_baro, func={"p" : self.get_3x3to6}
         ))
      dset(self,"p", depend_array(name='p', value=np.zeros((3,3),float),
            synchro=sync_baro, func={"p6" : self.get_6to3x3}
         ))

      if not p is None:
         self.p = p
      else:
         self.p = 0.0

      if not h0 is None:
         self.h0 = h0
      else:
         self.h0 = Cell()

      if not stressext is None:
         self.stressext = stressext
      else:
         self.stressext[:] = -1.0

   def bind(self, beads, nm, cell, forces, bias=None, prng=None, fixdof=None):
      """Binds beads, cell and forces to the barostat.

         This takes a beads object, a cell object and a forcefield object and
         makes them members of the barostat. It also then creates the objects that
         will hold the data needed in the barostat algorithms and the dependency
         network.

         Args:
         beads: The beads object from which the bead positions are taken.
         nm: The normal modes propagator object
         cell: The cell object from which the system box is taken.
         forces: The forcefield object from which the force and virial are
         taken.
         prng: The parent PRNG to bind the thermostat to
         fixdof: The number of blocked degrees of freedom.
         """

      super(BaroRGB, self).bind(beads, nm, cell, forces, bias, prng, fixdof)

      # obtain the thermostat mass from the given time constant (1/3 of what used for the corresponding NPT case)
      # note that the barostat temperature is nbeads times the physical T
      dset(self,"m", depend_array(name='m', value=np.atleast_1d(0.0),
                                 func=(lambda:np.asarray([self.tau**2*self.beads.natoms*Constants.kb*self.temp])),
                                 dependencies=[ dget(self,"tau"), dget(self,"temp") ] ))

      dset(self,"m6", depend_array(name='m6', value=np.zeros(6,float),
                                 func=(lambda:np.asarray([1,1,1,1,1,1])*self.m[0]),
                                 dependencies=[ dget(self,"m")] ))

      # overrides definition of pot to depend on the many things it depends on for anisotropic cell
      dset(self,"pot",
         depend_value(name='pot', func=self.get_pot,
            dependencies=[ dget(self.cell,"h"), dget(self.h0,"h"),
               dget(self.h0,"V"), dget(self.h0,"ih"), dget(self,"stressext") ]))

      # binds the thermostat to the piston degrees of freedom
      self.thermostat.bind(pm=[ self.p6, self.m6], prng=prng)

      dset(self,"kin",depend_value(name='kin',
            func=(lambda:0.5*np.trace(np.dot(self.p.T,self.p))/self.m[0]),
            dependencies= [dget(self,"p"), dget(self,"m")] ) )

      # the barostat energy must be computed from bits & pieces (overwrite the default)
      dset(self, "ebaro", depend_value(name='ebaro', func=self.get_ebaro,
                           dependencies=[ dget(self, "kin"), dget(self, "pot"),
                           dget(self.cell, "h"), dget(self, "temp"),
                           dget(self.thermostat,"ethermo")] ))

   def get_3x3to6(self):
      rp=np.zeros(6,float)
      rp[0]=self.p[0,0]; rp[1]=self.p[1,1]; rp[2]=self.p[2,2];
      rp[3]=self.p[0,1]; rp[4]=self.p[0,2]; rp[5]=self.p[1,2];
      return rp

   def get_6to3x3(self):
      rp=np.zeros((3,3),float)
      rp[0,0]=self.p6[0]; rp[1,1]=self.p6[1]; rp[2,2]=self.p6[2];
      rp[0,1]=self.p6[3]; rp[0,2]=self.p6[4]; rp[1,2]=self.p6[5];
      return rp

   def get_pot(self):
      """Calculates the elastic strain energy of the cell."""

      # NOTE: since there are nbeads replicas of the unit cell, the enthalpy contains a nbeads factor
      eps=np.dot(self.cell.h, self.h0.ih)
      eps=np.dot(eps.T, eps)
      eps=0.5*(eps - np.identity(3))

      return self.h0.V*np.trace(np.dot(self.stressext,eps))*self.beads.nbeads

   def get_ebaro(self):
      """Calculates the barostat conserved quantity."""

      lastterm=np.sum([(3-i)*np.log(self.cell.h[i][i]) for i in range(3)])
      lastterm = Constants.kb*self.temp*lastterm
      return self.thermostat.ethermo + self.kin + self.pot - lastterm

   def pstep(self):
      """Propagates the momenta for half a time step."""

      dthalf = self.dt*0.5
      dthalf2 = dthalf**2
      dthalf3 = dthalf**3/3.0

      hh0=np.dot(self.cell.h, self.h0.ih)
      pi_ext=np.dot(hh0, np.dot(self.stressext, hh0.T))*self.h0.V/self.cell.V
      L=np.diag([3,2,1])

      # This differs from the BZP thermostat in that it uses just one kT in the propagator.
      # This leads to an ensemble equaivalent to Martyna-Hughes-Tuckermann for both fixed and moving COM
      # Anyway, it is a small correction so whatever.

      self.p += dthalf*( self.cell.V* np.triu( self.stress - self.beads.nbeads*pi_ext ) +
                           Constants.kb*self.temp*L)

      fc = np.sum(depstrip(self.forces.f),0).reshape(self.beads.natoms,3)/self.beads.nbeads
      if self.bias != None: fc+= np.sum(depstrip(self.bias.f),0).reshape(self.beads.natoms,3)/self.beads.nbeads
      fcTonm = (fc/depstrip(self.beads.m3)[0].reshape(self.beads.natoms,3)).T
      pc = depstrip(self.beads.pc).reshape(self.beads.natoms,3)

      # I am not 100% sure, but these higher-order terms come from integrating the pressure virial term,
      # so they should need to be multiplied by nbeads to be consistent with the equations of motion in the PI context
      # again, these are tiny tiny terms so whatever.
      self.p += np.triu(dthalf2*np.dot(fcTonm,pc) + dthalf3*np.dot(fcTonm,fc)) * self.beads.nbeads

      self.beads.p += depstrip(self.forces.f)*dthalf
      if self.bias != None:  self.beads.p += depstrip(self.bias.f)*dthalf

   def qcstep(self):
      """Propagates the centroid position and momentum and the volume."""

      v = self.p/self.m[0]
      expq, expp = (matrix_exp(v*self.dt), matrix_exp(-v*self.dt))

      m = depstrip(self.beads.m)

      saveq=self.nm.qnm[0].copy()
      savep=self.nm.pnm[0].copy()
      for i in range(self.beads.natoms):
         self.nm.qnm[0,3*i:3*(i+1)] = np.dot(expq, self.nm.qnm[0,3*i:3*(i+1)])
         self.nm.qnm[0,3*i:3*(i+1)] += np.dot(np.dot(invert_ut3x3(v),(expq-expp)/(2.0)),depstrip(self.nm.pnm)[0,3*i:3*(i+1)]/m[i])
         self.nm.pnm[0,3*i:3*(i+1)] = np.dot(expp, self.nm.pnm[0,3*i:3*(i+1)])

      self.cell.h = np.dot(expq,self.cell.h)
