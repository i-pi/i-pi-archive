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

__all__ = ['Initializer', 'InitFile']

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


      for (k,v) in self.queue:
         ratoms = []   #r* means the new values from the initializer, i.e.
         rcell = None  # the values 'read' from the restart file.
         rbeads = Beads(0,0)

         info(" # inizializer parsing " + str(k) + " object.", verbosity.high)
         if k == "file"  or k == "file_v" or k == "file_p":
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

         if k == "cell":
            rcell = v

            if icell.V > 0.0:
               warning("Initialize from <cell> overwrites previous cell configuration.", verbosity.low)

            if rcell.V > 0.0:
               icell.h = rcell.h
            else:
               ValueError("Could not initialize the cell configuration from <initialize>.")


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
