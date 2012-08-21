"""Holds the class which computes important properties of the system, and
prepares them for output.

Classes:
   Properties: This is the class that holds all the algorithms to calculate
      the important properties that should be output.
   Trajectories: This class deals with outputting all position data in the
      appropriate format.
"""

__all__ = ['Properties', 'Trajectories', 'getkey']

import numpy as np
import math, random
from utils.depend import *
from utils.units import Constants, unit_to_internal, unit_to_user
from utils.mathtools import logsumlog, h2abc
from utils.io import *
from atoms import *
from cell import *
from ensembles import *
from forces import *

def getkey(pstring):
   """ Strips units and argument lists from a property/trajectory keyword. """

   pa=pstring.find('(');
   if pa<0: pa=len(pstring)
   pu=pstring.find('{');
   if pu<0: pu=len(pstring)
   return pstring[0:min(pa,pu)]


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
   """

   def __init__(self):
      """Initialises Properties."""

      self.property_dict = {}

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
      'density': Density of the system,
      'volume': Simulation box volume,
      'h': Cell vector matrix. Requires arguments x and v to give h[x,v],
      'potential': Potential energy estimator,
      'spring': The spring potential energy estimator,
      'kinetic_md': Classical kinetic energy estimator,
      'kinetic_cv': Quantum centroid virial kinetic energy estimator,
      'stress_md': The classical stress tensor estimator. Requires arguments
         x and v, to give stress[x,v],
      'stress_cv': The quantum centroid virial estimator of
         the stress tensor. Requires arguments x and v, to give stress[x,v],
      'pressure_md': Classical pressure estimator,
      'pressure_cv': Quantum centroid virial pressure estimator,
      'kstress_md': Classical kinetic stress tensor estimator. Requires
         arguments x and v, to give kstress[x,v],
      'kstress_cv': Quantum centroid virial kinetic stress tensor estimator.
         Requires arguments x and v, to give kstress[x,v],
      'virial_md': Classical virial tensor estimator. Requires arguments x and
         v, to give vir[x,v],
      'virial_cv': Quantum centroid virial virial tensor estimator. Requires
         arguments x and v, to give vir[x,v],
      'gle_ke': Kinetic energy for the additional momenta for the normal
         mode kinetic energy estimator. Takes one argument, which gives the
         mode that the kinetic energy is given for,
      'kin_yama': Quantum scaled coordinate kinetic energy estimator,
      'linlin': The scaled Fourier transform of the momentum distribution.
         Given by n(x) in Lin Lin et al., Phys. Rev. Lett. 105, 110602,
      'isotope_yama': Tcv(m/m') and log(R(m/m')) for isotope substitution
         calculations using the Yamamoto estimator,
      'isotope_thermo': Tcv(m/m') and log(R(m/m')) for isotope substitution
         calculations using the thermodynamic estimator.}.


      Args:
         simul: The Simulation object to be bound.
      """

      self.ensemble = simul.ensemble
      self.beads = simul.beads
      self.nm = simul.nm
      self.cell = simul.cell
      self.forces = simul.forces
      self.simul = simul


      self.property_dict["step"] = { "dimension" : "number", "func" : (lambda: (1 + self.simul.step)) }
      self.property_dict["time"] = { "dimension": "time", "func": (lambda: (1 + self.simul.step)*self.ensemble.dt) }
      self.property_dict["conserved"] = {"dimension" : "energy", "func" : self.get_econs }
      self.property_dict["temperature"] = { "dimension" : "temperature", "func" : self.get_temp }
      self.property_dict["density"] = {"dimension" : "density", "func": (lambda: self.beads.m.sum()/self.cell.V) }
      self.property_dict["volume"] = {"dimension": "volume", "func" :(lambda: self.cell.V) }
      self.property_dict["h"] = {"dimension" : "length", "func" :self.wrap_cell }

      self.property_dict["potential"] =  {"dimension" : "energy", "func": (lambda: self.forces.pot/self.beads.nbeads ) }
      self.property_dict["spring"] =     {"dimension" : "energy", "func": (lambda: self.beads.vpath*self.nm.omegan2) }
      self.property_dict["kinetic_md"] = {"dimension" : "energy", "func": (lambda: self.beads.kin/self.beads.nbeads) }
      self.property_dict["kinetic_cv"] = {"dimension" : "energy", "func": self.get_kincv }

      self.property_dict["stress_md"] = {"func" : self.get_stress }
      self.property_dict["pressure_md"] = {"func" : self.get_press }
      self.property_dict["kstress_md"] = {"func" : self.get_kstress }
      self.property_dict["virial_md"] = {"func" : self.get_vir }

      self.property_dict["stress_cv"] =   {"func" : self.get_stresscv      }
      self.property_dict["pressure_cv"] = {"func" : self.get_presscv       }
      self.property_dict["kstress_cv"] =  {"func" :  self.get_kstresscv    }
      self.property_dict["virial_cv"] =   {"func" : self.get_vircv         }

      self.property_dict["gle_ke"] =   {"func" : self.get_gleke             }

      self.property_dict["kin_yama"] = {"func" : self.get_kinyama           }

      self.property_dict["linlin"] =   {"func" : self.get_linlin            }
      self.property_dict["isotope_sc"] = {"func" : self.get_isotope_yama ,
        "help" :  "Scaled coordinates free energy perturbation scaled mass KE estimator. Prints everything which is needed to compute the kinetic energy for a isotope-substituted system. The 7 elements are: <h> <h^2> <T_CV> <T_CV^2> ln(<e^-h>) ln(|<T_CV e^-h>|) sign(<T_CV e^-h>). Mixed units, so outputs only in a.u." }

      self.property_dict["isotope_thermo"] = {"dimension" : "undefined", "func" : self.get_isotope_thermo, "size" : 7,
        "help" : "Thermodynamic free energy perturbation scaled mass KE estimator. Prints everything which is needed to compute the kinetic energy for a isotope-substituted system. The 7 elements are: <h> <h^2> <T_CV> <T_CV^2> ln(<e^-h>) ln(|<T_CV e^-h>|) sign(<T_CV e^-h>). Mixed units, so outputs only in a.u." }

      # dummy beads and forcefield objects so that we can use scaled and
      # displaced path estimators without changing the simulation bead
      # coordinates
      self.dbeads = simul.beads.copy()
      self.dforces = ForceBeads()
      self.dforces.bind(self.dbeads, self.simul.cell,  self.simul._forcemodel)

   def __getitem__(self, key):
      """Retrieves the item given by key.

      Note that if the key contains a string (arg1; arg2; ... )
      then it will add the appropriate arguments and value pairs
      to the calculation function of the property. Note the brackets and
      the semi-colon separators.

      Args:
         key: A string contained in property_dict.

      Returns:
         The property labelled by the keyword key.
      """

      args = []
      unit=""
      arglist=()
      unstart=len(key)
      argstart=unstart

      if '}' in key:
         # the property has a user-defined unit
         unstart = key.find('{')
         unstop = key.find('}', unstart)
         if unstop == -1:
            raise ValueError("Incorrect format in property units " + key)
         unit = key[unstart+1:unstop]
      if '(' in key:
         # If the property has additional arguments
         argstart = key.find('(')
         argstop = key.find(')', argstart)
         if argstop == -1:
            raise ValueError("Incorrect format in property arguments " + key)

         argstr = key[argstart:argstop+1]
         arglist = io_xml.read_tuple(argstr, delims="()", split=";", arg_type=str)


      key = key[0:min(unstart,argstart)] # strips the arguments from key name
      pkey = self.property_dict[key]

      if "dimension" in pkey and unit != "":
         return  unit_to_user(pkey["dimension"], unit, pkey["func"](*arglist))
      else:   return pkey["func"](*arglist)


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

   def get_kstress(self, x=0, v=0):
      """Calculates the classical kinetic stress tensor.

      Returns kstress[x,v].
      """

      return self.beads.kstress[x,v]

   def get_vir(self, x=0, v=0):
      """Calculates the classical virial tensor.

      Returns vir[x,v].
      """

      return self.forces.vir[x,v]

   def kstress_cv(self):
      """Calculates the quantum centroid virial kinetic stress tensor
      estimator.
      """

      kst = np.zeros((3,3),float)
      q = depstrip(self.beads.q)
      qc = depstrip(self.beads.qc)
      na3 = 3*self.beads.natoms
      for b in range(self.beads.nbeads):
         for i in range(3):
            for j in range(i,3):
               kst[i,j] += np.dot(q[b,i:na3:3] - qc[i:na3:3],
                  depstrip(self.forces.f[b])[j:na3:3])

      kst *= -1.0/float(self.beads.nbeads)
      for i in range(3):
         kst[i,i] += Constants.kb*self.ensemble.temp*(self.beads.natoms)
      return kst

   def get_kstresscv(self, x=0, v=0):
      """Calculates the quantum centroid virial kinetic stress tensor
      estimator.

      Returns kstress[x,v].
      """

      x = int(x)
      v = int(v)
      return self.kstress_cv()[x,v]

   def get_stresscv(self, x=0, v=0):
      """Calculates the quantum centroid virial stress tensor estimator.

      Returns stress[x,v].
      """

      x = int(x)
      v = int(v)
      stress = (self.forces.vir/float(self.beads.nbeads) + self.kstress_cv())/self.cell.V
      return stress[x,v]

   def get_presscv(self):
      """Calculates the quantum centroid virial pressure estimator."""

      return np.trace(self.forces.vir/float(self.beads.nbeads) + self.kstress_cv())/(3.0*self.cell.V)

   def get_vircv(self, x=0, v=0):
      """Calculates the quantum centroid virial tensor estimator.

      Returns vir[x,v].
      """

      return self.forces.vir[x,v]/float(self.beads.nbeads)

   def get_kincv(self, atom=""):
      """Calculates the quantum centroid virial kinetic energy estimator."""

      try:
         iatom=int(atom)
         latom=""
      except:
         iatom=-1
         latom=atom


      q=depstrip(self.beads.q)
      qc=depstrip(self.beads.qc)
      f=depstrip(self.forces.f)

      acv = 0.0
      for i in range(self.beads.natoms):
         if (atom != "" and iatom != i and latom != self.beads.names[i]):
            continue

         kcv = 0.0
         for b in range(self.beads.nbeads):
            kcv += np.dot(q[b,3*i:3*(i+1)] - qc[3*i:3*(i+1)], f[b,3*i:3*(i+1)])
         kcv *= -0.5/self.beads.nbeads
         kcv += 1.5*Constants.kb*self.ensemble.temp
         acv += kcv

      return acv

   def get_gleke(self, mode=0):
      """Calculates the kinetic energy of the nm-gle additional momenta.

      Args:
         mode: Gives the index of the normal mode thermostat we want the
            kinetic energy for.
      """

      gleke = 0.0
      s = depstrip(self.ensemble.thermostat.s)
      if len(s.shape) < 2:
         raise NameError("gle_ke called without a gle thermostat")
      elif len(s.shape) == 2:
         for alpha in range(len(s[0,:])):
            gleke += s[mode, alpha]**2/2.0
      else:
         for i in range(len(s[0,:,0])):
            for alpha in range(len(s[0,0,:])):
               gleke += s[mode, i, alpha]**2/2.0
      return gleke

   _DEFAULT_FINDIFF = 1e-5
   _DEFAULT_FDERROR = 1e-9
   _DEFAULT_MINFID = 1e-12
   def get_kinyama(self, fd_delta= - _DEFAULT_FINDIFF):
      """Calculates the quantum scaled coordinate kinetic energy estimator.

      Args:
         fd_delta: the relative finite difference in temperature to apply in computing
            finite-difference quantities. If it is negative, will be scaled down automatically
            to avoid discontinuities in the potential.
      Uses a finite difference method to calculate the kinetic energy estimator
      without requiring the forces as for the centroid virial estimator.
      """

      dbeta = abs(fd_delta)

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

         if (fd_delta < 0 and abs((vplus + vminus)/(v0*2) - 1.0) > _DEFAULT_FDERROR and dbeta > _DEFAULT_MINFID):
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

   def get_linlin(self, ux=0, uy=0, uz=0, atom='H'):
      """Gives the estimator for the momentum distribution, by opening the
      ring polymer path.

      Args:
         ux: The x component of the opening vector.
         uy: The y component of the opening vector.
         uz: The z component of the opening vector.
         atom: The atom type for which the path will be opened.
      """

      u = np.array([float(ux), float(uy), float(uz)])
      count = 0
      n = 0.0
      for at in range(self.beads.natoms):
         if self.beads.names[at] == atom:
            index = at
            # Keeps record of one of the atom indices so we can get
            # the mass later.

            count += 1
            self.dbeads.q[:] = self.beads.q[:]
            for bead in range(self.beads.nbeads):
               self.dbeads.q[bead,3*at:3*(at+1)] += self.opening(bead)*u
            n += math.exp(-(self.dforces.pot - self.forces.pot)/(self.ensemble.ntemp*Constants.kb*self.beads.nbeads))

      if count == 0:
         print "Warning, no atom with the given name found for lin-lin momentum distribution estimator."
         return 0.0
      else:
         n0 = math.exp(self.beads.m[index]*np.dot(u,u)*self.ensemble.temp*Constants.kb/(2*Constants.hbar**2))
         return n*n0/float(count)

   def wrap_cell(self, x=0, v=0):
      """Returns the the x-th component of the v-th cell vector."""

      return self.cell.h[x,v]

   def get_isotope_yama(self, alpha="1.0", atom=""):
      """Gives the components of the yamamoto scaled-mass KE estimator
      for a given atom index.

      Args:
         alpha: m'/m the mass ratio
         atom: the index of the atom to compute the isotope fractionation
            pair for, or a label

      Returns:
         a tuple from which one can reconstruct all that is needed to
         compute the SMKEE, and its statistical accuracy:
         (sum_deltah, sum_ke, log(sum(weights)), log(sum(weight*ke)),
            sign(sum(weight*ke)) )
      """

      try:
         iatom = int(atom)
         latom = ""
      except:
         iatom = -1
         latom = atom

      alpha = float(alpha)

      atcv = 0.0
      atcv2 = 0.0
      alogr = 0.0
      alogr2 = 0.0
      law = 0.0
      lawke = 0.0
      sawke = 1.0
      ni = 0

      # strips dependency control since we are not gonna change the true beads in what follows
      q = depstrip(self.beads.q)
      f = depstrip(self.forces.f)
      qc = depstrip(self.beads.qc)

      for i in range(self.beads.natoms):
         # selects only the atoms we care about
         if (atom != "" and iatom != i and latom != self.beads.names[i]):
            continue

         ni += 1

         # arranges coordinate-scaled beads in a auxiliary beads object
         self.dbeads.q[:] = q[:]
         for b in range(self.beads.nbeads):
            self.dbeads.q[b,3*i:3*(i+1)] = ( qc[3*i:3*(i+1)]+
                        np.sqrt(1.0/alpha)*(q[b,3*i:3*(i+1)]-qc[3*i:3*(i+1)])  )

         tcv = 0.0
         for b in range(self.beads.nbeads):
            tcv += np.dot( (self.dbeads.q[b,3*i:3*(i+1)]-self.dbeads.qc[3*i:3*(i+1)]),
                          self.dforces.f[b,3*i:3*(i+1)] )
         tcv *= -0.5/self.simul.nbeads
         tcv += 1.5*Constants.kb*self.simul.ensemble.temp

         logr = (self.dforces.pot-self.forces.pot)/(Constants.kb*self.simul.ensemble.temp*self.beads.nbeads)

         atcv += tcv
         atcv2 += tcv*tcv

         alogr += logr
         alogr2 += logr*logr;

         #accumulates log averages in a way which preserves accuracy
         if (ni==1):
            law = -logr
         else:
            (law, drop) = logsumlog( (law,1.0), (-logr,1.0))

         #here we need to take care of the sign of tcv, which might as well be negative... almost never but...
         if (ni==1):
            lawke = -logr + np.log(abs(tcv))
            sawke = np.sign(tcv);
         else:
            (lawke, sawke) = logsumlog( (lawke, sawke), (-logr+np.log(abs(tcv)), np.sign(tcv)) )

         print "CHECK", ni, logr, tcv, law, lawke
      if ni==0:
         raise ValueError("Couldn't find an atom which matched the argument of isotope_y")

      return (alogr, alogr2, atcv, atcv2, law, lawke, sawke)

   def get_isotope_thermo(self, alpha="1.0", atom=""):
      """Gives the components of the thermodynamic scaled-mass KE
      estimator for a given atom index.

      Args:
         alpha: m'/m the mass ratio
         atom: the index of the atom to compute the isotope fractionation
            pair for, or a label

      Returns:
         a tuple from which one can reconstruct all that is needed to
         compute the SMKEE:
         (sum_deltah, sum_ke, log(sum(weights)), log(sum(weight*ke)),
            sign(sum(weight*ke)) )
      """

      try:
         iatom = int(atom)
         latom = ""
      except:
         iatom = -1
         latom = atom

      alpha = float(alpha)

      atcv = 0.0
      alogr = 0.0
      atcv2 = 0.0
      alogr2 = 0.0
      law = 0.0
      lawke = 0.0
      sawke = 1.0
      ni = 0

      # strips dependency control since we are not gonna change the true beads in what follows
      q = depstrip(self.beads.q)
      f = depstrip(self.forces.f)
      qc = depstrip(self.beads.qc)

      for i in range(self.beads.natoms):
         # selects only the atoms we care about
         if (atom != "" and iatom != i and latom != self.beads.names[i]):
            continue

         ni += 1;

         spr = 0.0
         for b in range(1,self.beads.nbeads):
            for j in range(3*i,3*(i+1)):
               spr += (q[b,j]-q[b-1,j])**2
         for j in range(3*i,3*(i+1)):
            spr += (q[self.beads.nbeads-1,j]-q[0,j])**2

         spr *= 0.5*self.beads.m[i]*self.nm.omegan2

         # centroid virial contribution from atom i
         tcv = 0.0
         for b in range(self.beads.nbeads):
            tcv += np.dot( (q[b,3*i:3*(i+1)]-qc[3*i:3*(i+1)]), f[b,3*i:3*(i+1)])
         tcv *= -0.5/self.simul.nbeads
         tcv += 1.5*Constants.kb*self.simul.ensemble.temp

         logr = (alpha-1)*spr/(Constants.kb*self.simul.ensemble.temp*self.beads.nbeads)

         atcv += tcv
         atcv2 += tcv*tcv
         alogr += logr
         alogr2 += logr*logr

         #accumulates log averages in a way which preserves accuracy
         if (ni==1):
            law = -logr
         else:
            (law, drop) = logsumlog( (law,1.0), (-logr,1.0))

         #here we need to take care of the sign of tcv, which might as well be negative... almost never but...
         if (ni==1):
            lawke = -logr+np.log(abs(tcv))
            sawke=np.sign(tcv)
         else:
            (lawke, sawke) = logsumlog( (lawke, sawke), (-logr+np.log(abs(tcv)), np.sign(tcv)) )

      if ni==0:
         raise ValueError("Couldn't find an atom which matched the argument of isotope_y")

      return np.asarray([alogr, alogr2, atcv, atcv2, law, lawke, sawke])

class Trajectories(dobject):
   """A simple class to take care of output of trajectory data.

   Attributes:
      format: The file format for the output files.
      simul: The simulation object from which the position data will be
         obtained.
      fatom: A dummy beads object used so that individual replica trajectories
         can be output.
   """

   def __init__(self):
      """Initialises a Trajectories object.  """

      self.traj_dict = {}

   def bind(self, simul):
      """ Binds to a simulation object to fetch atomic and force data.

      Args:
         simul: The simulation object that will be managed by this Trajectories.
      """

      self.simul = simul
      self.fatom = simul.beads[0].copy()


      self.traj_dict["positions"] =  { "dimension" : "length", "func" : (lambda : 1.0*self.simul.beads.q) }
      self.traj_dict["velocities"] =  { "dimension" : "velocity", "func" : (lambda : self.simul.beads.p/self.simul.beads.m3) }
      self.traj_dict["forces"] =  { "dimension" : "force", "func" : (lambda : 1.0*self.simul.force.f) }
      self.traj_dict["kinetic_cv"] =  { "dimension" : "energy", "func" : self.get_akcv }
      self.traj_dict["kinetic_od"] =  { "dimension" : "energy", "func" : self.get_akcv_od }
      self.traj_dict["springs"] =  { "dimension" : "energy", "func" : self.get_aspr }
      self.traj_dict["r_gyration"] =  { "dimension" : "length", "func" : (lambda : 1.0*self.simul.beads.rg) }
      self.traj_dict["x_centroid"] =  { "dimension" : "length", "func" : (lambda : 1.0*self.simul.beads.qc)  }
      self.traj_dict["v_centroid"] =  { "dimension" : "length", "func" : (lambda : self.simul.beads.pc/self.simul.beads.m3[0])  }


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
      # helper arrays to make it more obvious what we are computing
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

   def get_aspr(self):
      """Calculates the contribution to the kinetic energy due to each degree
      of freedom.
      """

      rv = np.zeros(self.simul.beads.natoms*3)
      for b in range(1,self.simul.beads.nbeads):
         rv[:] += (self.simul.beads.q[b]-self.simul.beads.q[b-1])*(self.simul.beads.q[b]-self.simul.beads.q[b-1])*self.simul.beads.m3[b]
      rv[:] += (self.simul.beads.q[0]-self.simul.beads.q[self.simul.beads.nbeads-1])*(self.simul.beads.q[0]-self.simul.beads.q[self.simul.beads.nbeads-1])*self.simul.beads.m3[0]

      rv *= 0.5*self.simul.nm.omegan2
      return rv

   def __getitem__(self, key):
      """ Gets one of the trajectories. """

      args = []
      unit=""
      arglist=()
      unstart=len(key)
      argstart=unstart


      if '}' in key:
         # the property has a user-defined unit
         unstart = key.find('{')
         unstop = key.find('}', unstart)
         if unstop == -1:
            raise ValueError("Incorrect format in trajectory units " + key)
         unit = key[unstart+1:unstop]
      if '(' in key:
         # If the property has additional arguments
         argstart = key.find('(')
         argstop = key.find(')', argstart)
         if argstop == -1:
            raise ValueError("Incorrect format in trajectory arguments " + key)

         argstr = key[argstart:argstop+1]
         arglist = io_xml.read_tuple(argstr, delims="()", split=";", arg_type=str)


      key = key[0:min(unstart,argstart)] # strips the arguments from key name
      pkey = self.traj_dict[key]

      if "dimension" in pkey and unit != "":
         return  unit_to_user(pkey["dimension"], unit, 1.0) * pkey["func"](*arglist)
      else:   return pkey["func"](*arglist)


   def print_traj(self, what, stream, b=0, format="pdb"):
      """Prints out a frame of a trajectory for the specified quantity and bead.

      Args:
         what: A string specifying what to print.
         b: The bead index. Defaults to 0.
         stream: A reference to the stream on which data will be printed.
      """

      cq = self[what]
      if getkey(what) in [ "positions", "velocities", "forces" ] :
         self.fatom.q[:]= cq[b]
      else: self.fatom.q[:] = cq

      if format == "pdb":
         io_pdb.print_pdb(self.fatom, self.simul.cell, stream, title=("Traj: %s Step:  %10d  Bead:   %5d " % (what, self.simul.step+1, b) ) )
      elif format == "xyz":
         io_xyz.print_xyz(self.fatom, self.simul.cell, stream, title=("Traj: %s Step:  %10d  Bead:   %5d " % (what, self.simul.step+1, b) ) )
