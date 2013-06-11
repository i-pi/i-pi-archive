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
from utils.mathtools import logsumlog, h2abc_deg
from utils.io import *
from atoms import *
from cell import *
from ensembles import *
from forces import *
from utils.messages import verbosity, info

def getkey(pstring):
   """Strips units and argument lists from a property/trajectory keyword.

   Args:
      pstring: The string input by the user that specifies an output,
         which in general will specify units and argument lists.

   Returns: A string giving the keyword for the property, stripped of the
      argument lists and units key words.
   """

   pa = pstring.find('(')
   if pa < 0:
      pa = len(pstring)
   pu = pstring.find('{')
   if pu < 0:
      pu = len(pstring)
   return pstring[0:min(pa,pu)].strip()

def getall(pstring):
   """Returns the keyword, units and argument list separately.

   Args:
      pstring: The string input by the user that specifies an output,
         which in general will specify units and argument lists.

   Returns: A tuple giving the keyword for the property, and its units
      argument list and key word argument list.  
   """

   unit = ""
   arglist = ()
   kwarglist = {}
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
      for arg in arglist:
         # If a keyword argument is used
         equals = arg.find('=')
         if equals >= 0:
            kwarglist[arg[0:equals].strip()] = arg[equals+1:].strip()
            arglist = tuple(a for a in arglist if not a == arg)

   pstring = pstring[0:min(unstart,argstart)].strip() # strips the arguments from pstring name

   return (pstring, unit, arglist, kwarglist)

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
                      "help" : "The current simulation time step.",
                      'func': (lambda: (1 + self.simul.step))},
      "time": {       "dimension": "time",
                      "help": "The elapsed simulation time.",
                      'func': (lambda: (1 + self.simul.step)*self.ensemble.dt)},
      "conserved": {  "dimension": "energy",
                      "help": "The value of the conserved energy quantity per bead.",
                      'func': (lambda: self.ensemble.econs/float(self.beads.nbeads))},
      "temperature": {"dimension": "temperature",
                      "help": "The current physical temperature. Takes an argument 'atom', which can be either an atom label or index to specify which species to find the temperature of. If not specified, all atoms are used.",
                      'func': self.get_temp },
      "density": {    "dimension": "density",
                      "help": "The physical system density.",
                      'func': (lambda: self.beads.m.sum()/self.cell.V)},
      "volume": {     "dimension": "volume",
                      "help": "The volume of the cell box.",
                      'func': (lambda: self.cell.V) },
      "cell_h": {    "dimension" : "length",
                      "help": "Gives cell vector matrix. Returns the 6 components in the form [xx, yy, zz, xy, xz, yz].",
                      "size": 6,
                      "func": (lambda: self.flatten(self.cell.h))},
      "cell_abcABC": {"dimension" : "undefined",
                      "help": "Gives the lengths of the cell vectors and the angles between them in degrees as a list of the form [a, b, c, A, B, C], where A is the angle between the sides of length b and c in degrees, and B and C are defined similarly. Since there are a mixture of different units, these can only be output in atomic-units.",
                      "size": 6,
                      'func': (lambda: np.asarray(h2abc_deg(self.cell.h)))},
      "potential": {  "dimension" : "energy",
                      "help": "The physical system potential energy.",
                      'func': (lambda: self.forces.pot/self.beads.nbeads)},
      "spring": {     "dimension" : "energy",
                      "help": "The spring potential energy between the beads.",
                      'func': (lambda: self.beads.vpath*self.nm.omegan2)},
      "r_gyration": { "dimension" : "length",
                      "help" : "Gives the average radius of gyration of the selected atoms. Takes an argument 'atom', which can be either an atom label or index to specify which species to find the radius of gyration of. If not specified, all atoms are used.",
                      "func": self.get_rg},
      "kinetic_md":  {"dimension" : "energy",
                      "help": "The simulation kinetic energy. Takes an argument 'atom', which can be either an atom label or index to specify which species to find the temperature of. If not specified, all atoms are used.",
                      'func': self.get_kinmd},
      "kinetic_cv":  {"dimension" : "energy",
                      "help": "The physical system kinetic energy. Takes an argument 'atom', which can be either an atom label or index to specify which species to find the temperature of. If not specified, all atoms are used.",
                      'func': self.get_kincv},
      "kinetic_tens":{"dimension" : "energy",
                      "help" : "The physical system kinetic energy tensor. Returns the 6 components in the form [xx, yy, zz, xy, xz, yz]. Takes an argument 'atom', which can be either an atom label or index to specify which species to find the temperature of. If not specified, all atoms are used.",
                      "size" : 6,
                      "func" : self.get_ktens},
      "kinetic_ij":  {"dimension" : "energy",
                      "help" : "The physical system kinetic energy tensor, amongst atoms i and j. Returns the 6 components in the form [xx, yy, zz, xy, xz, yz]. Takes arguments 'i' and 'j', which give the indices of the appropriate two atoms.",
                      "size" : 6,
                      "func" : self.get_kij},
      "atom_x": {     "dimension" : "length",
                      "help": "Prints to properties the position (x,y,z) of a particle given its index. Takes arguments index and bead. If bead is not specified, refers to the centroid.",
                      "size" : 3,
                      "func" : (lambda atom="", bead="-1": self.get_atom_vec(self.beads.q, atom=atom, bead=bead))},
      "atom_v": {     "dimension" : "velocity",
                      "help": "Prints to properties the velocity (x,y,z) of a particle given its index. Takes arguments index and bead. If bead is not specified, refers to the centroid.",
                      "size" : 3,
                      "func" : (lambda atom="", bead="-1": self.get_atom_vec(self.beads.p/self.beads.m3, atom=atom, bead=bead))},
      "atom_p": {     "dimension" : "momentum",
                      "help": "Prints to properties the momentum (x,y,z) of a particle given its index. Takes arguments index and bead. If bead is not specified, refers to the centroid.",
                      "size" : 3,
                      "func" : (lambda atom="", bead="-1": self.get_atom_vec(self.beads.p, atom=atom, bead=bead))},
      "atom_f": {     "dimension" : "force",
                      "help": "Prints to properties the force (x,y,z) acting on a particle given its index. Takes arguments index and bead. If bead is not specified, refers to the centroid.",
                      "size" : 3,
                      "func" : (lambda atom="", bead="-1": self.get_atom_vec(self.forces.f, atom=atom, bead=bead))},
      "stress_md": {  "dimension": "pressure",
                      "size" : 6,
                      "help": "The simulation stress tensor. Returns the 6 components in the form [xx, yy, zz, xy, xz, yz].",
                      "func": (lambda: self.flatten((self.forces.vir + self.nm.kstress)/self.cell.V))},
      "pressure_md": {"dimension": "pressure",
                      "help": "The simulation pressure.",
                      "func": (lambda: np.trace((self.forces.vir + self.nm.kstress)/(3.0*self.cell.V)))},
      "kstress_md":  {"dimension": "pressure",
                      "size" : 6,
                      "help": "The simulation kinetic stress tensor. Returns the 6 components in the form [xx, yy, zz, xy, xz, yz].",
                      "func": (lambda: self.flatten(self.nm.kstress/self.cell.V))},
      "virial_md": {  "dimension": "pressure",
                      "size" : 6,
                      "help": "The simulation virial tensor. Returns the 6 components in the form [xx, yy, zz, xy, xz, yz].",
                      "func": (lambda: self.flatten(self.forces.vir/self.cell.V))},
      "stress_cv": {  "dimension": "pressure",
                      "size" : 6,
                      "help": "The physical system stress tensor.  Returns the 6 components in the form [xx, yy, zz, xy, xz, yz].",
                      "func": (lambda: self.flatten(self.forces.vir + self.kstress_cv())/(self.cell.V*self.beads.nbeads))},
      "pressure_cv": {"dimension": "pressure",
                      "help": "The physical pressure of the system.",
                      "func": (lambda: np.trace(self.forces.vir + self.kstress_cv())/(3.0*self.cell.V*self.beads.nbeads))},
      "kstress_cv":  {"dimension": "pressure",
                      "size" : 6,
                      "help": "The physical system kinetic stress tensor. Returns the 6 components in the form [xx, yy, zz, xy, xz, yz].",
                      "func": (lambda: self.flatten(self.kstress_cv()/(self.cell.V*self.beads.nbeads)))},
      "virial_cv": {  "dimension": "pressure",
                      "size" : 6,
                      "help": "The physical system virial stress tensor. Returns the 6 components in the form [xx, yy, zz, xy, xz, yz].",
                      "func": (lambda: self.flatten(self.forces.vir/(self.cell.V*self.beads.nbeads)))},
      "yamamoto": {   "help": "Gives the estimators required to calculate the Yamamoto finite difference approximation to the kinetic energy and constant volume heat capacity. Returns eps_v and eps_v', as defined in Takeshi M. Yamamoto, Journal of Chemical Physics, 104101, 123 (2005). As the two estimators have a different dimension, this can only be output in atomic units. Takes one argument, 'fd_delta', which gives the value of the finite difference parameter used. It defaults to " + str(-self._DEFAULT_FINDIFF) + ". If the value of 'fd_delta' is negative, then its magnitude will be reduced automatically by the code if the finite difference error becomes too large.",
                      'func': self.get_yama_estimators,
                      "size": 2},
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
      self.dforces.bind(self.dbeads, self.simul.cell,  self.simul.flist)

   def __getitem__(self, key):
      """Retrieves the item given by key.

      Note that if the key contains a string (arg1; arg2; ... )
      then it will pass the appropriate positional arguments to the
      calculation function of the property. Note the brackets and
      the semi-colon separators. If instead we have the syntax 
      (arg1=val1;arg2; ... ), then the keyword/value pair (arg1,val1)
      will be added to the keyword argument list. The appropriate key word
      arguments will then be passed to the calculation function instead.

      Similarly, if the key contains a string {unit}, then it will take
      the string 'unit' and use it to define the units that the property
      is output in.

      Args:
         key: A string contained in property_dict.

      Returns:
         The property labelled by the keyword key, along with its unit 
         keyword, and the argument lists for the function used to calculate
         the property specified by the keyword key.
      """

      (key, unit, arglist, kwarglist) = getall(key)
      pkey = self.property_dict[key]

      #pkey["func"](*arglist,**kwarglist) gives the value of the property 
      #in atomic units. unit_to_user() returns the value in the user 
      #specified units.
      if "dimension" in pkey and unit != "":
         return unit_to_user(pkey["dimension"], unit, pkey["func"](*arglist,**kwarglist))
      else:
         return pkey["func"](*arglist,**kwarglist)

   def get_atom_vec(self, prop_vec, atom="", bead="-1"):
      """Gives a vector for one atom.

      Args:
         prop_vec: An array from which to take the atomic vector from.
         atom: The index of the atom for which the vector will
            be output.
         bead: The index of the replica of the atom for which the
            vector will be output. If less than 0, then the centroid is used.
      """

      if atom == "":
         raise IndexError("Must specify the index for atom_vec property")
      atom = int(atom)
      bead = int(bead)
      if atom >= self.beads.natoms:
         raise IndexError("Cannot output atom_vec property as atom index %d is larger than the number of atoms" % atom)
      if bead >= self.beads.nbeads:
         raise IndexError("Cannot output atom_vec property as bead index %d is larger than the number of beads" % bead)

      if bead < 0:
         atom_vec = np.zeros(3)
         for b in range(self.beads.nbeads):
            atom_vec += prop_vec[b,3*atom:3*(atom+1)]
         return atom_vec/float(self.beads.nbeads)
      else:
         return prop_vec[bead,3*atom:3*(atom+1)]

   def flatten(self, vec_2D):
      """Takes a 3*3 upper-triangular array and returns it as a 1D array,
      of the form [xx, yy, zz, xy, xz, yz].
      """

      return np.array([vec_2D[0,0], vec_2D[1,1], vec_2D[2,2], vec_2D[0,1], vec_2D[0,2], vec_2D[1,2]])

   def get_rg(self, atom=""):
      """Calculates the radius of gyration of the ring polymers.

      Args:
         atom: If given, specifies the atom to give the gyration radius
            for. If not, the system average gyration radius is given.
      """

      try:
         #iatom gives the index of the atom to be studied
         iatom = int(atom)
         latom = ""
         if iatom >= self.beads.natoms:
            raise IndexError("Cannot output gyration radius as atom index %d is larger than the number of atoms" % iatom)
      except ValueError:
         #here 'atom' is a label rather than an index which is stored in latom
         iatom = -1
         latom = atom

      q = depstrip(self.beads.q)
      qc = depstrip(self.beads.qc)
      nat = self.beads.natoms
      nb = self.beads.nbeads
      rg_tot = 0.0
      ncount = 0
      for i in range(nat):
         if (atom != "" and iatom != i and latom != self.beads.names[i]):
            continue

         rg_at = 0.0
         for j in range(nb):
            dq = q[j,3*i:3*(i+1)] - qc[3*i:3*(i+1)]
            rg_at += np.dot(dq, dq)
         ncount += 1
         rg_tot += math.sqrt(rg_at/float(nb))

      if ncount == 0:
         raise IndexError("Couldn't find an atom which matched the argument of r_gyration")

      return rg_tot/float(ncount)

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
            if iatom >= self.beads.natoms:
               raise IndexError("Cannot output temperature as atom index %d is larger than the number of atoms" % iatom)
         except ValueError:
            #here 'atom' is a label rather than an index which is stored in latom
            iatom = -1
            latom = atom

         nat = 0
         for i in range(self.beads.natoms):
            if (iatom == i or latom == self.beads.names[i]): 
               nat += 1

         if nat == 0:
            raise IndexError("Couldn't find an atom which matched the argument of temperature")
         # "spreads" the COM removal correction evenly over all the atoms...
         kedof = self.get_kinmd(atom)/nat*(self.beads.natoms/(3.0*self.beads.natoms*self.beads.nbeads - mdof))

      return kedof/(0.5*Constants.kb)

   def kstress_cv(self):
      """Calculates the quantum centroid virial kinetic stress tensor
      estimator.

      Note that this is not divided by the volume or the number of beads.

      Returns:
         A 3*3 tensor with all the components of the tensor.
      """

      kst = np.zeros((3,3),float)
      q = depstrip(self.beads.q)
      qc = depstrip(self.beads.qc)
      pc = depstrip(self.beads.pc)
      m = depstrip(self.beads.m)
      fall = depstrip(self.forces.f)
      na3 = 3*self.beads.natoms

      for b in range(self.beads.nbeads):
         for i in range(3):
            for j in range(i,3):
               kst[i,j] -= np.dot(q[b,i:na3:3] - qc[i:na3:3],
                  fall[b,j:na3:3])

      # return the CV estimator MULTIPLIED BY NBEADS -- again for consistency with the virial, kstress_MD, etc...
      for i in range(3):
         kst[i,i] += self.beads.nbeads * ( np.dot(pc[i:na3:3],pc[i:na3:3]/m) )

      return kst

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
         if iatom >= self.beads.natoms:
            raise IndexError("Cannot output kinetic energy as atom index %d is larger than the number of atoms" % iatom)
      except ValueError:
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
         k = 3*i
         for b in range(self.beads.nbeads):
            kcv += (q[b,k] - qc[k])* f[b,k] + (q[b,k+1] - qc[k+1])* f[b,k+1] + (q[b,k+2] - qc[k+2])* f[b,k+2]
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
            if iatom >= self.beads.natoms:
               raise IndexError("Cannot output kinetic energy as atom index %d is larger than the number of atoms" % iatom)
         except ValueError:
            #here 'atom' is a label rather than an index which is stored in latom
            iatom = -1
            latom = atom

         pnm = depstrip(self.nm.pnm)
         dm3 = depstrip(self.nm.dynm3)
         kmd = 0.0
         for i in range(self.beads.natoms):
            if (atom != "" and iatom != i and latom != self.beads.names[i]):
               continue
            k = 3*i
            for b in range(self.beads.nbeads):
               kmd += (pnm[b,k]**2 + pnm[b,k+1]**2 + pnm[b,k+2]**2)/(2.0*dm3[b,k])
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
      if i >= self.beads.natoms:
         raise IndexError("Cannot output kinetic_ij as atom index %d is larger than the number of atoms" % i)
      if j >= self.beads.natoms:
         raise IndexError("Cannot output kinetic_ij as atom index %d is larger than the number of atoms" % j)
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
         if iatom >= self.beads.natoms:
            raise IndexError("Cannot output kinetic tensor as atom index %d is larger than the number of atoms" % iatom)
      except ValueError:
         #here 'atom' is a label rather than an index which is stored in latom
         iatom = -1
         latom = atom

      tkcv = np.zeros((6),float)
      for i in range(self.beads.natoms):
         if (atom != "" and iatom != i and latom != self.beads.names[i]):
            continue

         tkcv += self.get_kij(str(i), str(i))

      return tkcv

   def get_yama_estimators(self, fd_delta= - _DEFAULT_FINDIFF):
      """Calculates the quantum scaled coordinate kinetic energy estimator.

      Uses a finite difference method to calculate the estimators
      needed to calculate the energy and heat capacity of the system, as
      shown in Takeshi M. Yamamoto, Journal of Chemical Physics,
      104101, 123 (2005). Returns both eps_v and eps_v' as defined in
      the above article. Note that heat capacity is calculated as
      beta**2*kboltzmann*(<eps_v**2> - <eps_v>**2 - <eps_v'>), and the 
      energy of the system as <eps_v>.

      Args:
         fd_delta: the relative finite difference in temperature to apply in
         computing finite-difference quantities. If it is negative, will be
         scaled down automatically to avoid discontinuities in the potential.
      """

      dbeta = abs(float(fd_delta))
      beta = 1.0/(Constants.kb*self.ensemble.temp)

      qc = depstrip(self.beads.centroid.q)
      q = depstrip(self.beads.q)
      v0 = self.forces.pot/self.beads.nbeads
      while True:
         splus = math.sqrt(1.0 + dbeta)
         sminus = math.sqrt(1.0 - dbeta)

         for b in range(self.beads.nbeads):
            self.dbeads[b].q = qc*(1.0 - splus) + splus*q[b,:]
         vplus = self.dforces.pot/self.beads.nbeads

         for b in range(self.beads.nbeads):
            self.dbeads[b].q = qc*(1.0 - sminus) + sminus*q[b,:]
         vminus = self.dforces.pot/self.beads.nbeads

         if (fd_delta < 0 and abs((vplus + vminus)/(v0*2) - 1.0) > self._DEFAULT_FDERROR and dbeta > self._DEFAULT_MINFID):
            dbeta *= 0.5
            info("Reducing displacement in Yamamoto kinetic estimator", verbosity.low)
            continue
         else:
            eps = ((1.0 + dbeta)*vplus - (1.0 - dbeta)*vminus)/(2*dbeta)
            eps += 0.5*(3*self.beads.natoms)/beta

            eps_prime = ((1.0 + dbeta)*vplus + (1.0 - dbeta)*vminus - 2*v0)/(dbeta**2*beta)
            eps_prime -= 0.5*(3*self.beads.natoms)/beta**2

            break

      return np.asarray([eps, eps_prime])

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
         if iatom >= self.beads.natoms:
            raise IndexError("Cannot output scaled-mass kinetic energy estimator as atom index %d is larger than the number of atoms" % iatom)
      except ValueError:
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

      if ni == 0:
         raise IndexError("Couldn't find an atom which matched the argument of isotope_y")

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
         if iatom >= self.beads.natoms:
            raise IndexError("Cannot output scaled-mass kinetic energy estimator as atom index %d is larger than the number of atoms" % iatom)
      except ValueError:
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
         raise IndexError("Couldn't find an atom which matched the argument of isotope_y")

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
      # Note that here we want to return COPIES of the different arrays, so we make sure to make an operation in order not to return a reference.
      "positions": { "dimension" : "length",
                     "help": "Prints the coordinate trajectories.",
                     'func': (lambda : 1.0*self.simul.beads.q)},
      "velocities": {"dimension" : "velocity",
                     "help": "Prints the velocity trajectories.",
                     'func': (lambda : self.simul.beads.p/self.simul.beads.m3)},
      "momenta": {"dimension" : "momentum",
                     "help": "Prints the momentum trajectories.",
                     'func': (lambda : 1.0*self.simul.beads.p)},
      "forces": {    "dimension" : "force",
                     "help": "Prints the force trajectories.",
                     'func': (lambda : 1.0*self.simul.forces.f)},
      "extras": {    "dimension" : "",
                     "help": "Prints the extras trajectories.",
                     'func': (lambda : self.simul.forces.extras)},
      "kinetic_cv": {"dimension" : "energy",
                     "help": "Prints the kinetic energy for each bead, resolved into Cartesian components. Prints out a vector of the form [xx, yy, zz]",
                     'func': self.get_akcv},
      "kinetic_od": {"dimension" : "energy",
                     "help": "Prints the off diagonal elements of the kinetic stress tensor, for each bead. Prints out a vector of the form [xy, xz, yz]",
                     'func': self.get_akcv_od},
      "r_gyration": {"dimension" : "length",
                     "help": "Prints the radius of gyration for each atom, resolved into Cartesian components. Prints out a vector of the form [xx, yy, zz]",
                     'func': self.get_rg},
      "x_centroid": {"dimension" : "length",
                     "help": "Prints the centroid coordinates for each atom.",
                     'func': (lambda : 1.0*self.simul.beads.qc)},
      "v_centroid": {"dimension" : "velocity",
                     "help": "Prints the velocity centroid for each atom.",
                     'func': (lambda : self.simul.beads.pc/self.simul.beads.m3[0])},
      "p_centroid": {"dimension" : "momentum",
                     "help": "Prints the momentum centroid for each atom.",
                     'func': (lambda : 1.0*self.simul.beads.pc)},
      "f_centroid": {"dimension" : "force",
                     "help": "Prints the force centroid for each atom.",
                     'func': (lambda : np.sum(self.simul.forces.f,0)/float(self.simul.beads.nbeads))}
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
         rv[:,1] += dq[:,0]*f[:,2] + dq[:,2]*f[:,0]
         rv[:,2] += dq[:,1]*f[:,2] + dq[:,2]*f[:,1]
      rv *= 0.5
      rv *= -0.5/self.simul.beads.nbeads

      return rv.reshape(self.simul.beads.natoms*3)

   def get_rg(self):
      """Calculates the radius of gyration of the ring polymers.

      Computes separately the x, y, z contributions so that the actual
      gyration radius can be recovered as sqrt(rx^2+ry^2+rz^2).
      """

      q = depstrip(self.simul.beads.q)
      qc = depstrip(self.simul.beads.qc)
      nat = self.simul.beads.natoms
      nb = self.simul.beads.nbeads
      rg = np.zeros(3*nat)
      for i in range(nb):
         for j in range(nat):
            dq = q[i,3*j:3*(j+1)] - qc[3*j:3*(j+1)]
            rg[3*j:3*(j+1)] += dq*dq
      return np.sqrt(rg/float(nb))

   def __getitem__(self, key):
      """Retrieves the item given by key.

      Note that if the key contains a string (arg1; arg2; ... )
      then it will pass the appropriate positional arguments to the
      calculation function of the property. Note the brackets and
      the semi-colon separators. If instead we have the syntax 
      (arg1=val1;arg2; ... ), then the keyword/value pair (arg1,val1)
      will be added to the keyword argument list. The appropriate key word
      arguments will then be passed to the calculation function instead.

      Similarly, if the key contains a string {unit}, then it will take
      the string 'unit' and use it to define the units that the trajectory
      is output in.

      Args:
         key: A string contained in trajectory_dict.

      Returns:
         The trajectory labelled by the keyword key, along with its unit 
         keyword, and the argument lists for the function used to calculate
         the trajectory specified by the keyword key.
      """

      (key, unit, arglist, kwarglist) = getall(key)
      pkey = self.traj_dict[key]

      #pkey["func"](*arglist,**kwarglist) gives the value of the trajectory 
      #in atomic units. unit_to_user() returns the value in the user 
      #specified units.
      if "dimension" in pkey and unit != "":
         return  unit_to_user(pkey["dimension"], unit, 1.0) * pkey["func"](*arglist,**kwarglist)
      else:
         return pkey["func"](*arglist,**kwarglist)

   def print_traj(self, what, stream, b=0, format="pdb", cell_units="atomic_unit", flush=True):
      """Prints out a frame of a trajectory for the specified quantity and bead.

      Args:
         what: A string specifying what to print.
         b: The bead index. Defaults to 0.
         stream: A reference to the stream on which data will be printed.
         format: The output file format.
         cell_units: The units used to specify the cell parameters.
         flush: A boolean which specifies whether to flush the output buffer
            after each write to file or not.
      """

      cq = self[what]
      if getkey(what) in [ "extras" ] :
         stream.write(" #*EXTRAS*# Step:  %10d  Bead:  %5d  \n" % (self.simul.step+1, b) )
         stream.write(cq[b])
         stream.write("\n")
         return
      elif getkey(what) in [ "positions", "velocities", "forces" ] :
         self.fatom.q[:] = cq[b]
      else:
         self.fatom.q[:] = cq

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
