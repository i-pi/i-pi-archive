
from beads import Beads
from cell import Cell
from utils.io.io_xyz import read_xyz
from utils.io.io_pdb import read_pdb
from utils.depend import dobject

__all__ = ['Initializer', 'InitFile']

class InitFile(dobject):

   def __init__(self, filename="", format=""):

      self.filename=filename
      self.format=format


class Initializer(dobject):

   def __init__(self, nbeads=0, queue=None):

      if nbeads < 1: raise ValueError("Must specify number of beads when creating a Initializer object")
      self.nbeads = nbeads

      if queue is None:
         self.queue=[]
      else:
         self.queue=queue

   def init(self, simul):

      ibeads=Beads(0,0)
      icell=Cell()
      for (k,v) in self.queue:
         if k=="file" :
            # initialize from file

            rfile=open(v.filename,"r")
            ratoms=[]; rcell=None

            if (v.format=="xyz"):
               while True:
                  try:      myatoms = read_xyz(rfile)
                  except:   break
                  ratoms.append(myatoms)
            elif (v.format=="pdb"):
               while True:
                  try: myatoms, mycell = read_pdb(open(rfile,"r"))
                  except: break
                  ratoms.append(myatoms)
                  if rcell is None: rcell=mycell

            if not rcell is None:
               if icell.V > 0.0 :
                  print "WARNING: initialize from <file> overwrites previous cell configuration"
               icell=rcell

            if ibeads.nbeads > 0: print "WARNING: initialize from <file> overwrites previous path configuration"
            ibeads.resize(natoms=ratoms[0].natoms, nbeads=len(ratoms))
            ibeads.names=ratoms[0].names
            ibeads.m=ratoms[0].m
            for b in range(ibeads.nbeads):
               ibeads[b].q = ratoms[b].q
            # TODO implement other init modes!


      if ibeads.nbeads > 0:
         simul.beads.resize(ibeads.natoms, self.nbeads)
         # just start from bead number 0. must be fixed when path re-sampling is complete
         simul.beads.m=ibeads.m
         simul.beads.names=ibeads.names
         for b in range(self.nbeads):
            simul.beads.q[b] = ibeads[0].q; simul.beads.p[b] = ibeads[0].p
      elif self.beads.nbeads<1: raise ValueError("Could not initialize the path configuration, neither explicitly nor from <initialize>")

      # TODO work out a cleaner Cell object which can work for NPT, NVT etc
      #~ if icell.V > 0.0:
         #~ simul.cell=icell
      #~ elif simul.cell.V == 0.0 : raise ValueError("Could not initialize the cell configuration, neither explicitly nor from <initialize>")



