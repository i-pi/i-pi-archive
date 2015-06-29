#!/usr/bin/python
""" mergebeadspdb.py

Reads positions of individual beads from an i-PI run and 
assemles them in a pdb describing the ring polymer connectivity.

Assumes the input files are in pdb format names prefix.pos_*.pdb.

Syntax:
   mergebeadspdb.py prefix
"""

import numpy as np
import sys, glob
from ipi.utils.io import io_pdb
from ipi.engine.beads import Beads
from ipi.utils.depend import *
from ipi.utils.units import *

def main(prefix):

   ipos=[]
   for filename in sorted(glob.glob(prefix+".pos*")):
      ipos.append(open(filename,"r"))

   nbeads = len(ipos)   
   natoms = 0
   ifr = 0
   while True:
      try:
         for i in range(nbeads):
            pos, cell = io_pdb.read_pdb(ipos[i])
            if natoms == 0:
               natoms = pos.natoms
               beads = Beads(natoms,nbeads)
            beads[i].q = pos.q
      except EOFError: # finished reading files
         sys.exit(0)

      io_pdb.print_pdb_path(beads, cell)
      ifr+=1


if __name__ == '__main__':
   main(*sys.argv[1:])
