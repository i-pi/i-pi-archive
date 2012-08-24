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

from beads import Beads
from cell import Cell
from normalmodes import NormalModes
from ensembles import Ensemble

from utils.io.io_xyz import read_xyz
from utils.io.io_pdb import read_pdb
from utils.depend import dobject
from utils.units import Constants
from utils.nmtransform import nm_rescale
import numpy as np

__all__ = ['Initializer', 'InitFile']

class InitFile(dobject):

   def __init__(self, filename="", format=""):

      self.filename = filename
      self.format = format


class Initializer(dobject):

   def __init__(self, nbeads=0, queue=None):

      self.nbeads = nbeads

      if queue is None:
         self.queue = []
      else:
         self.queue = queue

   def init(self, simul):

      ibeads = simul.beads
      icell = simul.cell
      for (k,v) in self.queue:
         if k == "file" :
            # initialize from file

            rfile = open(v.filename,"r")
            ratoms = []
            rcell = None
            rbeads = Beads(0,0)
            if (v.format == "xyz"):
               while True:
                  try:
                     myatoms = read_xyz(rfile)
                  except:
                     break
                  ratoms.append(myatoms)
            elif (v.format == "pdb"):
               while True:
                  try:
                     myatoms, mycell = read_pdb(open(rfile,"r"))
                  except:
                     break
                  ratoms.append(myatoms)
                  if rcell is None:
                     rcell = mycell

            if not rcell is None:
               if icell.V > 0.0 :
                  print "WARNING: initialize from <file> overwrites previous cell configuration"
               icell = rcell

            rbeads.resize(natoms=ratoms[0].natoms, nbeads=len(ratoms))
            rbeads.names = ratoms[0].names
            rbeads.m = ratoms[0].m
            for b in range(rbeads.nbeads):
               rbeads[b].q = ratoms[b].q

            # scale rbeads up (or down) to self.nbeads!
            gbeads = Beads(rbeads.natoms,self.nbeads)
            res = nm_rescale(rbeads.nbeads,gbeads.nbeads)
            gbeads.q = res.b1tob2(rbeads.q)
            gbeads.m = rbeads.m
            gbeads.names = rbeads.names

            if ibeads.nbeads == self.nbeads:
               print "WARNING: initialize from <file> overwrites previous path configuration."
            else:
               ibeads.resize(rbeads.natoms,self.nbeads)

            if ibeads.natoms != gbeads.natoms:
               raise ValueError("Initialization tries to mix up structures with different atom numbers.")

            ibeads.q = gbeads.q
            ibeads.m = gbeads.m
            ibeads.names = gbeads.names

         if k == "beads":
            rbeads = v
            if rbeads.nbeads == self.nbeads:
               gbeads = rbeads
            else:
               gbeads = Beads(rbeads.natoms,self.nbeads)

               # scale rbeads up to self.nbeads!
               res = nm_rescale(rbeads.nbeads,gbeads.nbeads)
               gbeads.q = res.b1tob2(rbeads.q)
               gbeads.p = res.b1tob2(rbeads.p)

               gbeads.m = rbeads.m
               gbeads.names = rbeads.names

            if ibeads.nbeads > 0:
               print "WARNING: initialize from <beads> overwrites previous path configuration"
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
               for n in ibeads.names:
                  if n != "":
                     raise ValueError()
            except:
               ibeads.names = gbeads.names

         if k=="cell":
            pass
         # TODO work out a cleaner Cell object which can work for NPT, NVT etc
         #~ if icell.V > 0.0:
            #~ simul.cell=icell
         #~ elif simul.cell.V == 0.0 : raise ValueError("Could not initialize the cell configuration, neither explicitly nor from <initialize>")

         if k == "resample_v":
            if ibeads.natoms == 0:
               raise ValueError("Trying to resample velocities before having any structural information.")

            rtemp = v
            if rtemp < 0:
               rtemp = simul.ensemble.temp
            print "initializing at temperature", rtemp


            # pull together a mock initialization to get NM masses right 
            #without too much code duplication
            rbeads.resize(ibeads.natoms, ibeads.nbeads)
            rbeads.m[:] = ibeads.m
            rnm = NormalModes(mode=simul.nm.mode, freqs=simul.nm.nm_freqs)
            rens = Ensemble(dt=simul.ensemble.dt, temp=simul.ensemble.temp)
            rnm.bind(rbeads,rens)
            # then we exploit the sync magic to do a complicated initialization
            # in the NM representation
            # with (possibly) shifted-frequencies NM
            rnm.pnm = simul.prng.gvec((rbeads.nbeads,3*rbeads.natoms))*np.sqrt(rnm.dynm3)*np.sqrt(rbeads.nbeads*rtemp*Constants.kb)

            ibeads.p = rbeads.p

      print simul.beads.p
