"""Contains the classes that are used to initialize data in the simulation.

These classes can either be used to restart a simulation with some different
data or used to start a calculation. Any data given in these classes will
overwrite data given elsewhere.

Classes:
   Initializer: Holds the functions that are required to initialize objects in
      the code. Data can be initialized from a file, or according to a
      particular parameter. An example of the former would be initializing
      the configurations from a xyz file, an example of the latter would be
      initializing the velocities according to the physical temperature.
   InitFile: Simple class that allows initialization of data from a file.
      Just holds the information needed to open the file.
"""

import numpy as np

from beads import Beads
from cell import Cell
from normalmodes import NormalModes
from ensembles import Ensemble

from utils.io.io_xyz import read_xyz
from utils.io.io_pdb import read_pdb
from utils.io.io_xml import xml_parse_file
from utils.depend import dobject
from utils.units import Constants, unit_to_internal
from utils.nmtransform import nm_rescale
from utils.messages import verbosity, warning, info

__all__ = ['Initializer', 'InitFile', 'InitPositions', 'InitCell']

class InitFile(dobject):
   """Class that holds data about a particular input file.

   Attributes:
      filename: A string giving the name of the file.
      format: A string giving the extension of the file.
   """

   def __init__(self, filename="", format="", units="", cell_units=""):
      """Initializes InitFile.

      Args:
         filename: A string giving the name of the file.
         format: A string giving the extension of the file.
      """

      self.filename = filename
      self.format = format
      self.units = units
      self.cell_units = cell_units

class InitBase(dobject):
   """Base class for initializer objects.

   Attributes:
      value: A duck-typed stored value.
      mode: A string that determines how the value is to be interpreted.
   """

   def __init__(self, value=None, mode="", units=""):
      """Initializes InitFile.

      Args:
         filename: A string giving the name of the file.
         format: A string giving the extension of the file.
      """

      self.value = value
      self.mode = mode
      self.units = units

class InitVector(InitBase):
   """Vector initializer object.

   Attributes:
      value: A duck-typed stored value.
      mode: A string that determines how the value is to be interpreted.
   """

   def __init__(self, value=None, mode="", units="", index=-1, bead=-1):
      """Initializes InitFile.

      Args:
         filename: A string giving the name of the file.
         format: A string giving the extension of the file.
      """

      super(InitVector, self).__init__(value, mode, units)
      self.index = index
      self.bead = bead

class InitPositions(InitVector) : pass
class InitVelocities(InitVector) : pass
class InitMomenta(InitVector) : pass
class InitMasses(InitVector) : pass
class InitLabels(InitVector) : pass
class InitCell(InitBase): pass

def init_xyz(filename):

   rfile = open(filename,"r")
   ratoms = []
   while True:
   #while loop, so that more than one configuration can be given
   #so multiple beads can be initialized at once.
      try:
         myatoms = read_xyz(rfile)
      except:
         break
      ratoms.append(myatoms)
   return ratoms

def init_pdb(filename):

   rfile = open(filename,"r")
   ratoms = []
   while True:
   #while loop, so that more than one configuration can be given
   #so multiple beads can be initialized at once.
      try:
         myatoms, rcell  = read_xyz(rfile)
      except:
         break
      ratoms.append(myatoms)
   return ( ratoms, rcell ) # if multiple frames, the last cell is returned

def init_chk(filename):
   # reads configuration from a checkpoint file
   rfile = open(filename,"r")
   xmlchk = xml_parse_file(rfile) # Parses the file.

   from inputs.simulation import InputSimulation
   simchk = InputSimulation()
   simchk.parse(xmlchk.fields[0][1])
   rcell = simchk.cell.fetch()
   rbeads = simchk.beads.fetch()

   return (rbeads, rcell)

def init_beads(iif, nbeads):
   mode = iif.mode; value = iif.value
   if mode == "xyz" or mode == "pdb":
      if mode == "xyz": ratoms = init_xyz(value)
      if mode == "pdb": ratoms = init_pdb(value)[0]
      rbeads = Beads(ratoms[0].natoms,len(ratoms))
      for i in range(len(ratoms)): rbeads[i]=ratoms[i]
   elif mode == "chk":
      rbeads = init_chk(value)[0]
   elif mode == "manual":
      raise ValueError("Cannot initialize manually a whole beads object.")

   return rbeads

def init_vector(iif, nbeads, momenta=False):
   mode = iif.mode; value = iif.value
   if mode == "xyz" or mode == "pdb":
      rq = init_beads(iif, nbeads).q
   elif mode == "chk":
      if momenta: rq = init_beads(iif, nbeads).p
      else:       rp = init_beads(iif, nbeads).q
   elif mode == "manual":
      rq = value

   # determines the size of the input data
   if (rq.ndim > 1): # if the input has information on the number of beads
      nbeads = len(rq)
      natoms = len(rq[0])/3
   if mode == "manual":
      if iif.bead>=0: # if there is a bead specifier then we return a single bead slice
         nbeads = 1
      natoms = len(rq)/nbeads/3
      rq.shape = (nbeads,3*natoms)

   return rq

def set_vector(iif, dq, rq):
   (nbeads, natoms) = rq.shape; natoms/=3
   (dbeads, datoms) = dq.shape; datoms/=3

   # Check that indices make sense
   if iif.index < 0 and natoms!= datoms:
      raise ValueError("Initialization tries to mix up structures with different atom numbers.")
   if iif.index>= datoms:
      raise ValueError("Cannot initialize single atom as atom index %d is larger than the number of atoms" % iif.index)
   if iif.bead>= dbeads:
      raise ValueError("Cannot initialize single bead as bead index %d is larger than the number of beads" % iif.bead)

   if iif.bead < 0:   # we are initializing the path
      res = nm_rescale(nbeads,dbeads)  # path rescaler
      if nbeads != dbeads:
         warning(" # Initialize is rescaling from %5d beads to %5d beads" % (nbeads, dbeads), verbosity.low)
      if iif.index < 0:
         dq[:] = res.b1tob2(rq)
      else: # we are initializing a specific atom
         dq[:,3*iif.index:3*(iif.index+1)] = res.b1tob2(rq)
   else:  # we are initializing a specific bead
      if iif.index < 0:
         dq[iif.bead] = rq
      else:
         dq[iif.bead,3*iif.index:3*(iif.index+1)] = rq

class Initializer(dobject):
   """Class that deals with the initialization of data.

   This can either be used to initialize the atom positions and the cell data
   from a file, or to initialize them from a beads, atoms or cell object.

   Currently, we use a ring polymer contraction scheme to create a new beads
   object from one given in initialize if they have different numbers of beads,
   as described in the paper T. E. Markland and D. E. Manolopoulos, J. Chem.
   Phys. 129, 024105, (2008). If the new beads object has more beads than
   the beads object it was initialized from, we set the higher ring polymer
   normal modes to zero.

   Attributes:
      queue: A list of things to initialize. Each member of the list is a tuple
         of the form ('type', 'object'), where 'type' specifies what kind of
         initialization is being done, and 'object' gives the data to
         initialize it from.
   """

   def __init__(self, nbeads=0, queue=None):
      """Initializes Initializer.

      Arguments:
         nbeads: The number of beads that we need in the simulation. Not
            necessarily the same as the number of beads of the objects we are
            initializing the data from.
         queue: A list of things to initialize. Each member of the list is a
            tuple of the form ('type', 'object'), where 'type' specifies what
            kind of initialization is being done, and 'object' gives the data to
            initialize it from.
      """

      self.nbeads = nbeads

      if queue is None:
         self.queue = []
      else:
         self.queue = queue

   def init(self, simul):
      """Initializes the simulation.

      Takes a simulation object, and uses all the data in the initialization
      queue to fill up the data needed to run the simulation.

      Args:
         simul: A simulation object to be initialized.

      Raises:
         ValueError: Raised if there is a problem with the initialization,
            if something that should have been has not been, or if the objects
            that have been specified are not compatible with each other.
      """

      ibeads = simul.beads  #i* means the original values from the simulation
      icell = simul.cell    #object, i.e. the 'initial' values


      if simul.beads.nbeads == 0:
         fpos = fmom = fmass = flab = fcell = False   # we don't have an explicitly defined beads object
      else:
         fpos = fmom = fmass = flab = fcell = True
      for (k,v) in self.queue:
         info(" # Inizializer parsing " + str(k) + " object.", verbosity.high)

         if k == "cell":
            if fcell : warning("Overwriting previous cell parameters", verbosity.medium)
            if v.mode == "pdb":
               rh = init_pdb(v.value)[1].h
            elif v.mode == "chk":
               rh = init_chk(v.value)[0].cell.h
            else:
               rh = v.value.reshape((3,3))
            rh *= unit_to_internal("length",v.units,1.0)

            simul.cell.h = rh
            if simul.cell.V == 0.0:
               ValueError("Cell provided has zero volume")

            fcell = True
         elif k == "masses":
            if fmass : warning("Overwriting previous atomic masses", verbosity.medium)
            if v.mode == "manual":
               rm = v.value
            else:
               rm = init_beads(v, self.nbeads).m
            rm *= unit_to_internal("mass",v.units,1.0)

            if v.bead < 0:   # we are initializing the path
               if v.index < 0:
                  simul.beads.m = rm
               else: # we are initializing a specific atom
                  simul.beads.m[v.index:v.index+1] = rm
            else:
               raise ValueError("Cannot change the mass of a single bead")
            fmass = True

         elif k == "labels":
            if flab : warning("Overwriting previous atomic labels", verbosity.medium)
            if v.mode == "manual":
               rn = v.value
            else:
               rn = init_beads(v, self.nbeads).names

            if v.bead < 0:   # we are initializing the path
               if v.index < 0:
                  simul.beads.names = rn
               else: # we are initializing a specific atom
                  simul.beads.names[v.index:v.index+1] = rn
            else:
               raise ValueError("Cannot change the label of a single bead")
            flab = True

         elif k == "positions":
            if fpos : warning("Overwriting previous atomic positions", verbosity.medium)
            # read the atomic positions as a vector
            rq = init_vector(v, self.nbeads)
            rq *= unit_to_internal("length",v.units,1.0)
            (nbeads, natoms) = rq.shape;   natoms = natoms/3

            # check if we must initialize the simulation beads
            if simul.beads.nbeads == 0:
               if v.index >= 0: raise ValueError("Cannot initialize single atoms before the size of the system is known")
               simul.beads.resize(natoms,self.nbeads)

            set_vector(v, simul.beads.q, rq)
            fpos = True

         elif (k == "velocities" or k == "momenta") and v.mode == "thermal" :   # intercept here thermal initialization, so we don't need to check further down
            if fmom : warning("Overwriting previous atomic momenta", verbosity.medium)
            if simul.beads.natoms == 0:
               raise ValueError("Trying to resample velocities before having any structural information.")
            if not fmass:
               raise ValueError("Trying to resample velocities before having masses.")

            rtemp = v.value
            if rtemp <= 0:
               warning("Using the simulation temperature to resample velocities", verbosity.low)
               rtemp = simul.ensemble.temp
            else:
               warning(" # Resampling velocities at temperature %s" % rtemp, verbosity.low)

            # TODO -- Initialize a single atom!

            # pull together a mock initialization to get NM masses right
            #without too much code duplication
            rbeads = Beads(simul.beads.natoms, simul.beads.nbeads)
            rbeads.m[:] = simul.beads.m
            rnm = NormalModes(mode=simul.nm.mode, transform_method=simul.nm.transform_method, freqs=simul.nm.nm_freqs)
            rens = Ensemble(dt=simul.ensemble.dt, temp=simul.ensemble.temp)
            rnm.bind(rbeads,rens)
            # then we exploit the sync magic to do a complicated initialization
            # in the NM representation
            # with (possibly) shifted-frequencies NM
            rnm.pnm = simul.prng.gvec((rbeads.nbeads,3*rbeads.natoms))*np.sqrt(rnm.dynm3)*np.sqrt(rbeads.nbeads*rtemp*Constants.kb)

            simul.beads.p = rbeads.p
            fmom = True

         elif k == "momenta":
            if fmom : warning("Overwriting previous atomic momenta", verbosity.medium)
            # read the atomic momenta as a vector
            rp = init_vector(v, self.nbeads, momenta = True)
            rp *= unit_to_internal("momentum",v.units,1.0)
            (nbeads, natoms) = rp.shape;   natoms = natoms/3

            # checks if we must initialize the simulation beads
            if simul.beads.nbeads == 0:
               if v.index >=0 : raise ValueError("Cannot initialize single atoms before the size of the system is known")
               simul.beads.resize(natoms,self.nbeads)

            rp *= np.sqrt(self.nbeads/nbeads)
            set_vector(v, simul.beads.p, rp)
            fmom = True

         elif k == "velocities":
            if fmom : warning("Overwriting previous atomic momenta", verbosity.medium)
            # read the atomic velocities as a vector
            rv = init_vector(v, self.nbeads)
            rv *= unit_to_internal("velocity",v.units,1.0)
            (nbeads, natoms) = rv.shape;   natoms = natoms/3

            # checks if we must initialize the simulation beads
            if simul.beads.nbeads == 0 or not fmass:
               ValueError("Cannot initialize velocities before the masses of the atoms are known")
               simul.beads.resize(natoms,self.nbeads)

            warning(" # Initializing from velocities uses the previously defined masses -- not the masses inferred from the file -- to build momenta", verbosity.low)
            rv *= simul.beads.m3
            rv *= np.sqrt(self.nbeads/nbeads)
            set_vector(v, simul.beads.p, rv)
            fmom = True


         rcell=None
         ratoms=[]
         rbeads=Beads(0,0)
         if k == "file"  or k == "file_v" or k == "file_p":
            fmass = True
            # initialize from file (positions, velocities or momenta)
            # in this case 'v' is a InitFile instance.
            #! will do the first bit assuming we are reading positions,
            #! and then convert to the appropriate units further down

            rfile = open(v.filename,"r")
            if (v.format == "xyz"):
               while True:
               #while loop, so that more than one configuration can be given
               #so multiple beads can be initialized at once.
                  try:
                     myatoms = read_xyz(rfile)
                  except:
                     break
                  ratoms.append(myatoms)

            elif (v.format == "pdb"):
               while True:
               #while loop, so that more than one configuration can be given
               #so multiple beads can be initialized at once.
                  try:
                     myatoms, mycell = read_pdb(rfile)
                  except:
                     break
                  ratoms.append(myatoms)
                  if k == "file" and rcell is None:
                     mycell.h *= unit_to_internal("length",v.units,1.0)
                     rcell = mycell

            elif (v.format == "chk" or v.format == "checkpoint"):
               # reads configuration from a checkpoint file
               rfile = open(v.filename,"r")
               xmlchk = xml_parse_file(rfile) # Parses the file.

               if k == "file_v":
                  warning(" Reading from checkpoint actually initializes momenta, not velocities. Make sure this is what you want. ",
                          verbosity.low)

               from inputs.simulation import InputSimulation
               simchk = InputSimulation()
               simchk.parse(xmlchk.fields[0][1])
               if k == "file":
                  rcell = simchk.cell.fetch()
               rbeads = simchk.beads.fetch()

            if not rcell is None:
               warning(" Initialize from <file> overwrites previous cell configuration. ", verbosity.low)
               icell.h = rcell.h

            if not (v.format == "chk" or v.format == "checkpoint"):
               # assembles the list of atomic configurations into a beads object
               rbeads.resize(natoms=ratoms[0].natoms, nbeads=len(ratoms))
               rbeads.names = ratoms[0].names
               rbeads.m = ratoms[0].m
               if k == "file":
                  for b in range(rbeads.nbeads):
                     rbeads[b].q = ratoms[b].q * unit_to_internal("length",v.units,1.0)
               elif k == "file_p":
                  for b in range(rbeads.nbeads):
                     rbeads[b].p = ratoms[b].q * unit_to_internal("momentum",v.units,1.0)
               elif k == "file_v":
                  for b in range(rbeads.nbeads):
                     rbeads[b].p = ratoms[b].q * rbeads.m3 * unit_to_internal("velocity",v.units,1.0)


            # scale rbeads up (or down) to self.nbeads!
            gbeads = Beads(rbeads.natoms,self.nbeads)
            res = nm_rescale(rbeads.nbeads,gbeads.nbeads)
            if rbeads.nbeads != gbeads.nbeads:
               warning(" # Initialize is rescaling from %5d beads to %5d beads" % (rbeads.nbeads, self.nbeads),
                     verbosity.low)

            gbeads.q = res.b1tob2(rbeads.q)
            gbeads.p = res.b1tob2(rbeads.p) * np.sqrt(gbeads.nbeads/rbeads.nbeads)

            ### CAUTION! THIS MAY BE WRONG WHEN (DE)CONTRACTING THE RING POLYMER. SHOULD CHECK CAREFULLY!
            #TODO: is this checked? then the comment should go
            gbeads.m = rbeads.m
            gbeads.names = rbeads.names

            if ibeads.nbeads == self.nbeads:
               warning("Initialize from <file> overwrites previous path configuration.", verbosity.low)
            else:
               ibeads.resize(rbeads.natoms,self.nbeads)

            if ibeads.natoms != gbeads.natoms:
               raise ValueError("Initialization tries to mix up structures with different atom numbers.")

            if k == "file":
               ibeads.q = gbeads.q
               ibeads.m = gbeads.m
               ibeads.names = gbeads.names
            else:  # chk files always have the momenta, but we don't touch them unless required.
               ibeads.p = gbeads.p

         if k == "beads":
            rbeads = v
            info("names " + str(v.names), verbosity.high)
            info("vq" + str(v.m), verbosity.high)
            info("uq" + str(np.linalg.norm(v.q)), verbosity.high)

            if rbeads.nbeads == self.nbeads:
               gbeads = rbeads
            else:
               warning(" # Initialize is rescaling from %5d beads to %5d beads" % (rbeads.nbeads, self.nbeads),
                     verbosity.low)
               gbeads = Beads(rbeads.natoms,self.nbeads)

               # scale rbeads up to self.nbeads!
               res = nm_rescale(rbeads.nbeads,gbeads.nbeads)
               gbeads.q = res.b1tob2(rbeads.q)
               gbeads.p = res.b1tob2(rbeads.p) * np.sqrt(gbeads.nbeads/rbeads.nbeads)   ### CAUTION! THIS MAY BE WRONG WHEN (DE)CONTRACTING THE RING POLYMER. SHOULD CHECK CAREFULLY!

               gbeads.m = rbeads.m
               gbeads.names = rbeads.names

            if ibeads.nbeads > 0:
               warning("Initialize from <beads> overwrites previous path configuration.", verbosity.low)
            else:
               ibeads.resize(rbeads.natoms,self.nbeads)

            if ibeads.natoms != gbeads.natoms:
               raise ValueError("Initialization tries to mix up structures with different atom numbers.")

            # consider a vectors of zeros as a "ignore field" statement

            if np.linalg.norm(gbeads.q) > 0.0:
               ibeads.q = gbeads.q
            if np.linalg.norm(gbeads.p) > 0.0:
               ibeads.p = gbeads.p
            if np.linalg.norm(gbeads.m) > 0.0:
               ibeads.m = gbeads.m
            try:
               for n in gbeads.names:
                  if n != "":
                     raise ValueError()
            except:
               ibeads.names = gbeads.names




         if k == "resample_v":
            if ibeads.natoms == 0:
               raise ValueError("Trying to resample velocities before having any structural information.")

            rtemp = v
            if rtemp < 0:
               rtemp = simul.ensemble.temp

            # pull together a mock initialization to get NM masses right
            #without too much code duplication
            rbeads.resize(ibeads.natoms, ibeads.nbeads)
            rbeads.m[:] = ibeads.m
            rnm = NormalModes(mode=simul.nm.mode, transform_method=simul.nm.transform_method, freqs=simul.nm.nm_freqs)
            rens = Ensemble(dt=simul.ensemble.dt, temp=simul.ensemble.temp)
            rnm.bind(rbeads,rens)
            # then we exploit the sync magic to do a complicated initialization
            # in the NM representation
            # with (possibly) shifted-frequencies NM
            rnm.pnm = simul.prng.gvec((rbeads.nbeads,3*rbeads.natoms))*np.sqrt(rnm.dynm3)*np.sqrt(rbeads.nbeads*rtemp*Constants.kb)

            ibeads.p = rbeads.p

      if ibeads.natoms == 0:
         raise ValueError("Could not initialize the path configuration, neither explicitly nor from <initialize>")
      if icell.V == 0.0 :
         raise ValueError("Could not initialize the cell configuration, neither explicitly nor from <initialize>")
