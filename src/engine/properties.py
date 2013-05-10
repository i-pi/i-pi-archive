"""Holds the class which computes important properties of the system, and
prepares them for output.

Classes:
   Properties: This is the class that holds all the algorithms to calculate
      the important properties that should be output.
   Trajectories: This class deals with outputting all position data in the
      appropriate format.

Functions:
   getkey: This function strips the units and argument list specification
      from a string specifying an output parameter.
   getall: This function gives the keyword, units and argument list
      specification from a string specifying an output parameter.
   help_latex: This returns a string that can be used in the manual to
      specify the different available outputs.
"""

__all__ = ['Properties', 'Trajectories', 'getkey', 'getall', 'help_latex']

import os
import numpy as np
import math
from utils.depend import *
from utils.units import Constants, unit_to_internal, unit_to_user
from utils.mathtools import logsumlog, h2abc
from utils.io import *
from atoms import *
from cell import *
from ensembles import *
from forces import *

def getkey(pstring):
   """Strips units and argument lists from a property/trajectory keyword.

   Args:
      pstring: The string input by the user that specifies an output,
         which in general will specify units and argument lists.

   Returns: A string giving the keyword for the property.
   """

   pa = pstring.find('(')
   if pa < 0:
      pa = len(pstring)
   pu = pstring.find('{')
   if pu < 0:
      pu = len(pstring)
   return pstring[0:min(pa,pu)]

def getall(pstring):
   """Returns the keyword, units and argument list separately.

   Args:
      pstring: The string input by the user that specifies an output,
         which in general will specify units and argument lists.

   Returns: A tuple giving the keyword for the property, and its units
      and argument list.
   """

   unit = ""
   arglist = ()
   unstart = len(pstring)
   argstart = unstart

   if '}' in pstring:
      # the property has a user-defined unit
      unstart = pstring.find('{')
      unstop = pstring.find('}', unstart)
      if unstop == -1:
         raise ValueError("Incorrect format in units specification " + pstring)
      unit = pstring[unstart+1:unstop]
   if '(' in pstring:
      # If the property has additional arguments
      argstart = pstring.find('(')
      argstop = pstring.find(')', argstart)
      if argstop == -1:
         raise ValueError("Incorrect format in argument list " + pstring)

      argstr = pstring[argstart:argstop+1]
      arglist = io_xml.read_tuple(argstr, delims="()", split=";", arg_type=str)

   pstring = pstring[0:min(unstart,argstart)] # strips the arguments from pstring name

   return (pstring, unit, arglist)

def help_latex(idict, ref=False):
   """Function to generate a LaTeX formatted file.

   Args:
      idict: Either property_dict or traj_dict, to be used to
         generate the help file.
      ref: A boolean giving whether the latex file produced will be a
         stand-alone document, or will be intended as a section of a larger
         document with cross-references between the different sections.

   Returns:
      A LaTeX formatted string.
   """

   rstr = ""
   if not ref:
      #assumes that it is a stand-alone document, so must have document
      #options.
      rstr += "\\documentclass[12pt,fleqn]{report}"
      rstr += "\n\\begin{document}\n"
      rstr += "The following are the different allowable ouputs:\n"
   rstr += "\\begin{itemize}\n"

   for out in idict:
      rstr += "\\item {\\bf " + out + "}: " + idict[out]['help'] + "\n"
      try:
         if idict[out]['dimension'] != "undefined":
            #doesn't print out dimension if not necessary.
            dimstr = "\n {\\bf DIMENSION}: " + idict[out]['dimension'] + '\n'
            rstr += dimstr
      except KeyError:
         pass
      try:
         sizestr = "\n{\\bf SIZE}: " + str(idict[out]['size']) + '\n'
         rstr += sizestr
      except KeyError:
         pass

   rstr += "\\end{itemize}\n"

   if not ref:
      #ends the created document if it is not part of a larger document
      rstr += "\\end{document}"

   # Some escape characters are necessary for the proper latex formatting
   rstr = rstr.replace('_', '\\_')
   rstr = rstr.replace('\\\\_', '\\_')
   rstr = rstr.replace('...', '\\ldots ')
   rstr = rstr.replace('<', '$<$')
   rstr = rstr.replace('>', '$>$')
   rstr = rstr.replace('|', '$|$')

   return rstr


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
      nm: A normal modes object giving the normal mode representation.
      cell: A cell object giving the system box.
      forces: A forcefield object giving the force calculator for each
         replica of the system.
      property_dict: A dictionary containing all the properties that can be
         output.
   """

   _DEFAULT_FINDIFF = 1e-5
   _DEFAULT_FDERROR = 1e-9
   _DEFAULT_MINFID = 1e-12

   def __init__(self):
      """Initialises Properties."""

      self.property_dict = {
      "step": {       "dimension" : "number",
                      "help" : "The current time step of the simulation.",
                      'func': (lambda: (1 + self.simul.step))},
      "time": {       "dimension": "time",
                      "help": "The elapsed simulation time.",
                      'func': (lambda: (1 + self.simul.step)*self.ensemble.dt)},
      "conserved": {  "dimension": "energy",
                      "help": "The value of the conserved energy quantity per bead.",
                      'func': self.get_econs},
      "temperature": {"dimension": "temperature",
                      "help": "The current physical temperature.",
                      'func': self.get_temp },
      "density": {    "dimension": "density",
                      "help": "The physical density of the system.",
                      'func': (lambda: self.beads.m.sum()/self.cell.V)},
      "volume": {     "dimension": "volume",
                      "help": "The volume of the cell box.",
                      'func': (lambda: self.cell.V) },
      "cell_h": {     "dimension" : "length",
                      "help": "Gives one of the cell parameters. Takes arguments 'x' and 'v', which gives h[x,v]. By default gives h[0,0].",
                      'func': self.wrap_cell},
      "cell_h6": {    "dimension" : "length",
                      "help": "Gives all the non-zero cell parameters, in the order h_xx, h_yy, h_xz, h_xy, h_xz, h_yz.",
                      "size": 6,
                      'func': self.full_cell},
      "cell_abcABC": {"dimension" : "undefined",
                      "help": "Gives the lengths of the cell vectors and the angles between them in degrees as a list. Since there are a mixture of different units, these can only be output in atomic-units.",
                      "size": 6,
                      'func': self.cell_abcABC},
      "potential": {  "dimension" : "energy",
                      "help": "The potential energy of the system.",
                      'func': (lambda: self.forces.pot/self.beads.nbeads)},
      "spring": {     "dimension" : "energy",
                      "help": "The spring potential energy between the beads.",
                      'func': (lambda: self.beads.vpath*self.nm.omegan2)},
      "kinetic_md":  {"dimension" : "energy",
                      "help": "The classical kinetic energy of the simulation.",
                      'func': self.get_kinmd},
      "kinetic_cv":  {"dimension" : "energy",
                      "help": "The physical kinetic energy of the system.",
                      'func': self.get_kincv},
      "kinetic_tens":{"dimension" : "energy",
                      "help" : "The T_xx T_yy T_zz T_xy T_xz T_yz components of the kinetic energy tensor (c-v estimator).",
                      "size" : 6,
                      "func" : self.get_ktens},
      "kinetic_ij":  {"dimension" : "energy",
                      "help" : "The many-body T_xx T_yy T_zz T_xy T_xz T_yz components of the kinetic energy tensor, amongst atoms i and j (c-v estimator).",
                      "size" : 6,
                      "func" : self.get_kij},
      "atom_x": {     "dimension" : "length",
                      "help": "Prints to properties the position (x,y,z) of a particle given its index. Takes arguments index and bead. If bead is not specified, refers to the centroid.",
                      "size" : 3,
                      'func': self.get_atomx},
      "atom_v": {     "dimension" : "velocity",
                      "help": "Prints to properties the velocity (x,y,z) of a particle given its index. Takes arguments index and bead. If bead is not specified, refers to the centroid.",
                      "size" : 3,
                      'func': self.get_atomv},
      "stress_md": {  "dimension": "pressure",
                      "help": "The classical stress tensor of the simulation. Takes arguments 'x' and 'v', which gives stress[x,v]. By default gives stress[0,0].",
                      'func': self.get_stress},
      "pressure_md": {"dimension": "pressure",
                      "help": "The classical pressure of the simulation.",
                      'func': self.get_press},
      "kstress_md":  {"dimension": "pressure",
                      "help": "The classical kinetic stress tensor of the simulation. Takes arguments 'x' and 'v', which gives kstress[x,v]. By default gives kstress[0,0].",
                      'func': self.get_kstress},
      "virial_md": {  "dimension": "pressure",
                      "help": "The classical virial tensor of the simulation. Takes arguments 'x' and 'v', which gives virial[x,v]. By default gives virial[0,0].",
                      'func': self.get_vir},
      "stress_cv": {  "dimension": "pressure",
                      "help": "The physical stress tensor of the system. Takes arguments 'x' and 'v', which gives stress[x,v]. By default gives stress[0,0].",
                      'func': self.get_stresscv},
      "pressure_cv": {"dimension": "pressure",
                      "help": "The physical pressure of the system.",
                      'func': self.get_presscv},
      "kstress_cv":  {"dimension": "pressure",
                      "help": "The physical kinetic stress tensor of the system. Takes arguments 'x' and 'v', which gives kstress[x,v]. By default gives kstress[0,0].",
                      'func': self.get_kstresscv},
      "virial_cv": {  "dimension": "pressure",
                      "help": "The physical virial stress tensor of the system. Takes arguments 'x' and 'v', which gives virial[x,v]. By default gives virial[0,0].",
                      'func': self.get_vircv},
#Tensor versions of all of these as well?
      "gle_ke": {     "dimension": "energy",
                      "help": "Gives the kinetic energy associated with the additional degrees of freedom used in the GLE thermostat. Takes an argument 'mode' which gives the degree of freedom that is looked at, and defaults to 0.",
                      'func': self.get_gleke},
#This is currently horribly messy, and should probably be removed.
      "kin_yama": {   "dimension": "energy",
                      "help": "Gives the Yamamoto kinetic energy estimator. Takes one argument, 'fd_delta', which gives the value of the finite difference parameter used. It defaults to " + str(-self._DEFAULT_FINDIFF) + ".",
                      'func': self.get_kinyama},
      "isotope_scfep":  {"dimension": "undefined",
                      "size": 7,
                      'func': self.get_isotope_yama,
                      "help" :  "Scaled coordinates free energy perturbation scaled mass KE estimator. Prints everything which is needed to compute the kinetic energy for a isotope-substituted system. The 7 elements are: <h> <h**2> <T_CV> <T_CV**2> ln(sum(e**(-h))) ln(|sum(T_CV e**(-h))|) sign(sum(T_CV e**(-h))). Mixed units, so outputs only in a.u. Takes two arguments, 'alpha' and 'atom', which give the scaled mass parameter and the atom of interest respectively, and default to '1.0' and ''. The 'atom' argument can either be the label of a particular kind of atom, or an index of a specific atom." },
      "isotope_tdfep":  {"dimension" : "undefined",
                          "size" : 7,
                          'func': self.get_isotope_thermo,
                          "help" : "Thermodynamic free energy perturbation scaled mass KE estimator. Prints everything which is needed to compute the kinetic energy for a isotope-substituted system. The 7 elements are: <h> <h**2> <T_CV> <T_CV**2> ln(sum(e**(-h))) ln(|sum(T_CV e**(-h))|) sign(sum(T_CV e**(-h))). Mixed units, so outputs only in a.u. Takes two arguments, 'alpha' and 'atom', which give the scaled mass parameter and the atom of interest respectively, and default to '1.0' and ''. The 'atom' argument can either be the label of a particular kind of atom, or an index of a specific atom." }
      }

   def bind(self, simul):
      """Binds the necessary objects from the simulation to calculate the
      required properties.

      Args:
         simul: The Simulation object to be bound.
      """

      self.ensemble = simul.ensemble
      self.beads = simul.beads
      self.nm = simul.nm
      self.cell = simul.cell
      self.forces = simul.forces
      self.simul = simul
      # dummy beads and forcefield objects so that we can use scaled and
      # displaced path estimators without changing the simulation bead
      # coordinates
      self.dbeads = simul.beads.copy()
      self.dforces = Forces()
      self.dforces.bind(self.dbeads, self.simul.cell,  self.simul.flist, self.simul.soft_exit, self.simul.verb)

   def __getitem__(self, key):
      """Retrieves the item given by key.

      Note that if the key contains a string (arg1; arg2; ... )
      then it will add the appropriate arguments and value pairs
      to the calculation function of the property. Note the brackets and
      the semi-colon separators.

      Similarly, if the key contains a string {unit}, then it will take
      the string 'unit' and use it to define the units that the property
      is output in.

      Args:
         key: A string contained in property_dict.

      Returns:
         The property labelled by the keyword key.
      """

      (key, unit, arglist) = getall(key)
      pkey = self.property_dict[key]

      #pkey["func"](*arglist) gives the value of the property in atomic units
      #unit_to_user returns the value in the user specified units.
      if "dimension" in pkey and unit != "":
         return  unit_to_user(pkey["dimension"], unit, pkey["func"](*arglist))
      else:
         return pkey["func"](*arglist)

   def get_atomx(self, atom="", bead="-1"):
      """Gives the position vector of one atom.

      Args:
         atom: The index of the atom for which the position vector will
            be output.
         bead: The index of the replica of the atom for which the position
            vector will be output. If less than 0, then the centroid is used.
      """

      if atom == "":
         raise ValueError("Must specify the index for atom_x property")
      atom = int(atom)
      bead = int(bead)
      if bead < 0:
         return self.beads.centroid[atom].q
      else:
         return self.beads[bead][atom].q

   def get_atomv(self, atom="", bead="-1"):
      """Gives the velocity vector of one atom.

      Args:
         atom: The index of the atom for which the velocity vector will
            be output.
         bead: The index of the replica of the atom for which the velocity
            vector will be output. If less than 0, then the centroid is used.
      """

      if atom == "":
         raise ValueError("Must specify the index for atom_x property")
      atom = int(atom)
      bead = int(bead)
      if bead < 0:
         return self.beads.centroid[atom].p/ self.beads.m[atom]
      else:
         return self.beads[bead][atom].p/ self.beads.m[atom]

   def get_temp(self, atom=""):
      """Calculates the MD kinetic temperature.

      Note that in the case that the centre of mass constraint there will be
      3 fewer degrees of freedom than without, so this has to be taken into
      account when calculating the kinetic temperature.

      Args:
         atom: If given, specifies the atom to give the temperature
            for. If not, then the simulation temperature.
      """

      if self.ensemble.fixcom:
         mdof = 3
      else:
         mdof = 0


      if atom == "":
         # use the KE computed in the NM representation in order to avoid problems when mass scaling is used
         kedof = self.get_kinmd()/(3*self.beads.natoms*self.beads.nbeads - mdof)
      else:
         try:
            #iatom gives the index of the atom to be studied
            iatom = int(atom)
            latom = ""
         except:
            #here 'atom' is a label rather than an index which is stored in latom
            iatom = -1
            latom = atom

         nat = 0
         for i in range(self.beads.natoms):
            if (iatom == i or latom == self.beads.names[i]): nat+=1

         # "spreads" the COM removal correction evenly over all the atoms...
         kedof = self.get_kinmd(atom)/nat*(self.beads.natoms/(3.0*self.beads.natoms*self.beads.nbeads - mdof))

      return kedof/(0.5*Constants.kb)

   def get_econs(self):
      """Calculates the conserved quantity estimator per bead."""

      return self.ensemble.econs/(self.beads.nbeads)

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

      x = int(x)
      v = int(v)
      return self.beads.kstress[x,v]/self.cell.V

   def get_vir(self, x=0, v=0):
      """Calculates the classical virial tensor.

      Returns vir[x,v].
      """

      x = int(x)
      v = int(v)
      #TODO Maybe this should not be divided by V?
      return self.forces.vir[x,v]/self.cell.V

   def kstress_cv(self):
      """Calculates the quantum centroid virial kinetic stress tensor
      estimator.

      Note, in line with the beads.kstress tensor, this does not divide
      by the volume, and so this must be done before the results are
      output.
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
      return self.kstress_cv()[x,v]/self.cell.V

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

      x = int(x)
      v = int(v)
      return self.forces.vir[x,v]/(self.beads.nbeads*self.cell.V)

   def get_kincv(self, atom=""):
      """Calculates the quantum centroid virial kinetic energy estimator.

      Args:
         atom: If given, specifies the atom to give the kinetic energy
            for. If not, the system kinetic energy is given.
      """

      try:
         #iatom gives the index of the atom to be studied
         iatom = int(atom)
         latom = ""
      except:
         #here 'atom' is a label rather than an index which is stored in latom
         iatom = -1
         latom = atom


      q = depstrip(self.beads.q)
      qc = depstrip(self.beads.qc)
      f = depstrip(self.forces.f)

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

   def get_kinmd(self, atom=""):
      """Calculates the classical kinetic energy of the simulation (p^2/2m)

      Args:
         atom: If given, specifies the atom to give the kinetic energy
            for. If not, the simulation kinetic energy is given.
      """

      if atom == "":
         return self.nm.kin/self.beads.nbeads
      else:
         try:
            #iatom gives the index of the atom to be studied
            iatom = int(atom)
            latom = ""
         except:
            #here 'atom' is a label rather than an index which is stored in latom
            iatom = -1
            latom = atom

         pnm = depstrip(self.nm.pnm)
         dm3 = depstrip(self.nm.dynm3)
         kmd = 0.0
         for i in range(self.beads.natoms):
            if (atom != "" and iatom != i and latom != self.beads.names[i]):
               continue
            for b in range(self.beads.nbeads):
               kmd += (pnm[b,3*i]**2 + pnm[b,3*i+1]**2 + pnm[b,3*i+2]**2)/(2.0*dm3[b,3*i])
         return kmd/self.beads.nbeads

   def get_kij(self, ni="0", nj="0"):
      """Calculates the quantum centroid virial kinetic energy
      TENSOR estimator for two possibly different atom indices.

      Args:
         ni: The index of atom i.
         nj: The index of atom j.

      Returns:
         The contribution to the kinetic energy tensor estimator from
         the interactions between atom i and atom j.
      """

      i = int(ni)
      j = int(nj)
      mi = self.beads.m[i]
      mj = self.beads.m[j]
      ai = 3*i
      aj = 3*j

      q = depstrip(self.beads.q)
      qc = depstrip(self.beads.qc)
      f = depstrip(self.forces.f)

      # I implement this for the most general case. In practice T_ij = <p_i p_j>/(2sqrt(m_i m_j))
      kcv = np.zeros((6),float)
      for b in range(self.beads.nbeads):
         kcv[0] += mi*(q[b,ai] - qc[ai])    *f[b,aj]   + mj*(q[b,aj] - qc[aj])    *f[b,ai]       #Txx
         kcv[1] += mi*(q[b,ai+1] - qc[ai+1])*f[b,aj+1] + mj*(q[b,aj+1] - qc[aj+1])*f[b,ai+1]     #Tyy
         kcv[2] += mi*(q[b,ai+2] - qc[ai+2])*f[b,aj+2] + mj*(q[b,aj+2] - qc[aj+2])*f[b,ai+2]     #Tzz
         kcv[3] += mi*(q[b,ai] - qc[ai])*    f[b,aj+1] + mj*(q[b,aj+1] - qc[aj+1])*f[b,ai]       #Txy
         kcv[4] += mi*(q[b,ai] - qc[ai])*    f[b,aj+2] + mj*(q[b,aj+2] - qc[aj+2])*f[b,ai]       #Txz
         kcv[5] += mi*(q[b,ai+1] - qc[ai+1])*f[b,aj+2] + mj*(q[b,aj+2] - qc[aj+2])*f[b,ai+1]     #Tyz

      kcv *= -0.5/(self.beads.nbeads*2*np.sqrt(mi*mj))
      if i == j:
         kcv[0:3] += 0.5*Constants.kb*self.ensemble.temp

      return kcv


   def get_ktens(self, atom=""):
      """Calculates the quantum centroid virial kinetic energy
      TENSOR estimator.

      Args:
         atom: The index of the atom for which the kinetic energy tensor
            is to be output, or the index of the type of atoms for which
            it should be output.
      """

      try:
         #iatom gives the index of the atom to be studied
         iatom = int(atom)
         latom = ""
      except:
         #here 'atom' is a label rather than an index which is stored in latom
         iatom = -1
         latom = atom


      tkcv = np.zeros((6),float)
      for i in range(self.beads.natoms):
         if (atom != "" and iatom != i and latom != self.beads.names[i]):
            continue

         tkcv += self.get_kij(str(i), str(i))

      return tkcv

   def get_gleke(self, mode=0):
      """Calculates the kinetic energy of the nm-gle additional momenta.

      Args:
         mode: Gives the index of the normal mode thermostat we want the
            kinetic energy for.
      """

      mode = int(mode)
      gleke = 0.0
      try:
         s = depstrip(self.ensemble.thermostat.s)
      except AttributeError:
         return 0.0

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

   def get_kinyama(self, fd_delta= - _DEFAULT_FINDIFF):
      """Calculates the quantum scaled coordinate kinetic energy estimator.

      Uses a finite difference method to calculate the kinetic energy estimator
      without requiring the forces as for the centroid virial estimator.

      Args:
         fd_delta: the relative finite difference in temperature to apply in
         computing finite-difference quantities. If it is negative, will be
         scaled down automatically to avoid discontinuities in the potential.
      """

      dbeta = abs(float(fd_delta))

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

         if (fd_delta < 0 and abs((vplus + vminus)/(v0*2) - 1.0) > self._DEFAULT_FDERROR and dbeta > self._DEFAULT_MINFID):
            dbeta *= 0.5
            print "Reducing displacement in Yamamoto kinetic estimator"
            continue
         else:
            break

      return kyama

   def wrap_cell(self, x=0, v=0):
      """Returns the the x-th component of the v-th cell vector."""

      x = int(x)
      v = int(v)
      return self.cell.h[x,v]

   def full_cell(self):
      """Returns a six-component vector giving all the non-zero components of
      the cell vector matrix.
      """

      h = depstrip(self.cell.h)
      return np.array([h[0,0], h[1,1], h[2,2], h[0,1], h[0,2], h[1,2]])

   def cell_abcABC(self):
      """Returns the cell parameters as the length of the principle cell
      vectors and the angles between them."""

      h = depstrip(self.cell.h)
      (a, b, c, alpha, beta, gamma) = h2abc(h)
      return np.array([a, b, c, alpha*180/math.pi, beta*180/math.pi, gamma*180/math.pi])

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
         #iatom gives the index of the atom to be studied
         iatom = int(atom)
         latom = ""
      except:
         #here 'atom' is a label rather than an index which is stored in latom
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
         tcv *= -0.5/self.beads.nbeads
         tcv += 1.5*Constants.kb*self.simul.ensemble.temp

         logr = (self.dforces.pot-self.forces.pot)/(Constants.kb*self.simul.ensemble.temp*self.beads.nbeads)

         atcv += tcv
         atcv2 += tcv*tcv

         alogr += logr
         alogr2 += logr*logr;

         #accumulates log averages in a way which preserves accuracy
         if (ni == 1):
            law = -logr
         else:
            (law, drop) = logsumlog( (law,1.0), (-logr,1.0))

         #here we need to take care of the sign of tcv, which might as well be
         #negative... almost never but...
         if (ni == 1):
            lawke = -logr + np.log(abs(tcv))
            sawke = np.sign(tcv);
         else:
            (lawke, sawke) = logsumlog( (lawke, sawke), (-logr+np.log(abs(tcv)), np.sign(tcv)) )

         print "CHECK", ni, logr, tcv, law, lawke
      if ni == 0:
         raise ValueError("Couldn't find an atom which matched the argument of isotope_y")

      return np.asarray([alogr/ni, alogr2/ni, atcv/ni, atcv2/ni, law, lawke, sawke])

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
         #iatom gives the index of the atom to be studied
         iatom = int(atom)
         latom = ""
      except:
         #here 'atom' is a label rather than an index which is stored in latom
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

         ni += 1

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
         tcv *= -0.5/self.beads.nbeads
         tcv += 1.5*Constants.kb*self.simul.ensemble.temp

         logr = (alpha-1)*spr/(Constants.kb*self.simul.ensemble.temp*self.beads.nbeads)


         atcv += tcv
         atcv2 += tcv*tcv
         alogr += logr
         alogr2 += logr*logr

         #accumulates log averages in a way which preserves accuracy
         if (ni == 1):
            law = -logr
         else:
            (law, drop) = logsumlog( (law,1.0), (-logr,1.0))

         #here we need to take care of the sign of tcv, which might as well be
         #negative... almost never but...
         if (ni == 1):
            lawke = -logr + np.log(abs(tcv))
            sawke = np.sign(tcv)
         else:
            (lawke, sawke) = logsumlog( (lawke, sawke), (-logr+np.log(abs(tcv)), np.sign(tcv)) )

      if ni == 0:
         raise ValueError("Couldn't find an atom which matched the argument of isotope_y")

      return np.asarray([alogr/ni, alogr2/ni, atcv/ni, atcv2/ni, law, lawke, sawke])


class Trajectories(dobject):
   """A simple class to take care of output of trajectory data.

   Attributes:
      simul: The simulation object from which the position data will be
         obtained.
      fatom: A dummy beads object used so that individual replica trajectories
         can be output.
      traj_dict: A dictionary containing all the trajectories that can be
         output.
   """

   def __init__(self):
      """Initialises a Trajectories object."""

      self.traj_dict = {
      "positions": { "dimension" : "length",
                     "help": "Prints the coordinate trajectories.",
                     'func': (lambda : 1.0*self.simul.beads.q)},
      "velocities": {"dimension" : "velocity",
                     "help": "Prints the velocity trajectories.",
                     'func': (lambda : self.simul.beads.p/self.simul.beads.m3)},
      "forces": {    "dimension" : "force",
                     "help": "Prints the force trajectories.",
                     'func': (lambda : 1.0*self.simul.forces.f)},
      "extras": {    "dimension" : "",
                     "help": "Prints the extras trajectories.",
                     'func': (lambda : self.simul.forces.extras)},
      "kinetic_cv": {"dimension" : "energy",
                     "help": "Prints the kinetic energy for each bead, resolved into Cartesian components.",
                     'func': self.get_akcv},
      "kinetic_od": {"dimension" : "energy",
                     "help": "Prints the off diagonal elements of the kinetic stress tensor, for each bead.",
                     'func': self.get_akcv_od},
      "springs": {   "dimension" : "energy",
                     "help": "Prints the spring potential for each atom, resolved into Cartesian components.",
                     'func': self.get_aspr},
#This may be deprecated.
      "r_gyration": {"dimension" : "length",
                     "help": "Prints the radius of gyration for each atom.",
                     'func': (lambda : 1.0*self.simul.beads.rg)},
#r_gyration might be more suitable as a property.
      "x_centroid": {"dimension" : "length",
                     "help": "Prints the centroid coordinates for each atom.",
                     'func': (lambda : 1.0*self.simul.beads.qc)},
      "v_centroid": {"dimension" : "length",
                     "help": "Prints the velocity centroid for each atom.",
                     'func': (lambda : self.simul.beads.pc/self.simul.beads.m3[0])}
      }


   def bind(self, simul):
      """ Binds to a simulation object to fetch atomic and force data.

      Args:
         simul: The simulation object that will be managed by this Trajectories.
      """

      self.simul = simul
      self.fatom = simul.beads[0].copy()

   def get_akcv(self):
      """Calculates the contribution to the kinetic energy due to each degree
      of freedom.
      """

      rv = np.zeros(self.simul.beads.natoms*3)
      for b in range(self.simul.beads.nbeads):
         rv[:] += (self.simul.beads.q[b]-self.simul.beads.qc)*self.simul.forces.f[b]
      rv *= -0.5/self.simul.beads.nbeads
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
         rv[:,0] += dq[:,0]*f[:,1] + dq[:,1]*f[:,0]
         rv[:,1] += dq[:,1]*f[:,2] + dq[:,2]*f[:,1]
         rv[:,2] += dq[:,0]*f[:,2] + dq[:,2]*f[:,0]
      rv *= 0.5
      rv *= -0.5/self.simul.beads.nbeads

      return rv.reshape(self.simul.beads.natoms*3)

   def get_aspr(self):
      """Calculates the contribution to the kinetic energy due to each degree
      of freedom.
      """
      #TODO What the hell does this do?

      rv = np.zeros(self.simul.beads.natoms*3)
      for b in range(1,self.simul.beads.nbeads):
         rv[:] += (self.simul.beads.q[b]-self.simul.beads.q[b-1])*(self.simul.beads.q[b]-self.simul.beads.q[b-1])*self.simul.beads.m3[b]
      rv[:] += (self.simul.beads.q[0]-self.simul.beads.q[self.simul.beads.nbeads-1])*(self.simul.beads.q[0]-self.simul.beads.q[self.simul.beads.nbeads-1])*self.simul.beads.m3[0]

      rv *= 0.5*self.simul.nm.omegan2
      return rv

   def __getitem__(self, key):
      """ Gets one of the trajectories. """

      (key, unit, arglist) = getall(key)
      pkey = self.traj_dict[key]

      if "dimension" in pkey and unit != "":
         return  unit_to_user(pkey["dimension"], unit, 1.0) * pkey["func"](*arglist)
      else:
         return pkey["func"](*arglist)

   def print_traj(self, what, stream, b=0, format="pdb", cell_units="atomic_unit", flush=True):
      """Prints out a frame of a trajectory for the specified quantity and bead.

      Args:
         what: A string specifying what to print.
         b: The bead index. Defaults to 0.
         stream: A reference to the stream on which data will be printed.
      """

      cq = self[what]
      if getkey(what) in [ "extras" ] :
         stream.write(" #*EXTRAS*# Step:  %10d  Bead:  %5d  \n" % (self.simul.step+1, b) )
         stream.write(cq[b])
         stream.write("\n")
         return
      elif getkey(what) in [ "positions", "velocities", "forces" ] :
         self.fatom.q[:] = cq[b]
      else: self.fatom.q[:] = cq

      fcell = Cell()
      fcell.h = self.simul.cell.h*unit_to_user("length", cell_units, 1.0)

      if format == "pdb":
         io_pdb.print_pdb(self.fatom, fcell, stream, title=("Traj: %s Step:  %10d  Bead:   %5d " % (what, self.simul.step+1, b) ) )
      elif format == "xyz":
         io_xyz.print_xyz(self.fatom, fcell, stream, title=("Traj: %s Step:  %10d  Bead:   %5d " % (what, self.simul.step+1, b) ) )
      elif format == "bin":
         io_binary.print_bin(self.fatom, fcell, stream, title=("Traj: %s Step:  %10d  Bead:   %5d " % (what, self.simul.step+1, b) ) )
      if flush :
         stream.flush()
         os.fsync(stream)
