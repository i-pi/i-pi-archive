
from beads import Beads
from cell import Cell
from normalmodes import NormalModes
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

      self.nbeads = nbeads

      if queue is None:
         self.queue=[]
      else:
         self.queue=queue

   def init(self, simul):

      ibeads=simul.beads
      icell=simul.cell
      for (k,v) in self.queue:
         if k=="file" :
            # initialize from file

            rfile=open(v.filename,"r")
            ratoms=[]; rcell=None
            rbeads=Beads(0,0)
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

            rbeads.resize(natoms=ratoms[0].natoms, nbeads=len(ratoms))
            rbeads.names=ratoms[0].names
            rbeads.m=ratoms[0].m
            for b in range(rbeads.nbeads):
               rbeads[b].q = ratoms[b].q

            # TODO scale rbeads up to self.nbeads!
            gbeads=Beads(rbeads.natoms,self.nbeads)
            for b in range(self.nbeads): gbeads[b].q = rbeads[0].q
            gbeads.m=rbeads.m; gbeads.names=rbeads.names

            if ibeads.nbeads == self.nbeads: print "WARNING: initialize from <file> overwrites previous path configuration."
            else: ibeads.resize(rbeads.natoms,self.nbeads)

            if ibeads.natoms != gbeads.natoms: raise ValueError("Initialization tries to mix up structures with different atom numbers.")

            ibeads.q=gbeads.q
            ibeads.m=gbeads.m; ibeads.names=gbeads.names

         if k=="beads":
            rbeads=v
            if rbeads.nbeads == self.nbeads: gbeads=rbeads
            else:
               # TODO scale rbeads up to self.nbeads!
               gbeads=Beads(rbeads.natoms,self.nbeads)
               for b in range(self.nbeads):
                  gbeads[b].q = rbeads[0].q; gbeads[b].p = rbeads[0].p
               gbeads.m=rbeads.m; gbeads.names=rbeads.names

            if ibeads.nbeads > 0: print "WARNING: initialize from <beads> overwrites previous path configuration"
            else: ibeads.resize(rbeads.natoms,self.nbeads)

            if ibeads.natoms != gbeads.natoms: raise ValueError("Initialization tries to mix up structures with different atom numbers.")

            # consider a vectors of zeros as a "ignore field" statement
            if np.linalg.norm(gbeads.q) > 0.0: ibeads.q = gbeads.q
            if np.linalg.norm(gbeads.p) > 0.0: ibeads.p = gbeads.p
            if np.linalg.norm(gbeads.m) > 0.0: ibeads.m = gbeads.m
            try:
               for n in ibeads.names:
                  if n!="": raise ValueError()
            except: ibeads.names=gbeads.names

         if k=="cell": pass
         # TODO work out a cleaner Cell object which can work for NPT, NVT etc
         #~ if icell.V > 0.0:
            #~ simul.cell=icell
         #~ elif simul.cell.V == 0.0 : raise ValueError("Could not initialize the cell configuration, neither explicitly nor from <initialize>")

         if k=="resample":
            inm=NormalModes()
            pass




