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
from ipi.utils.io.backends.io_xyz import read_xyz
from ipi.utils.io.backends.io_pdb import print_pdb
from ipi.engine.beads import Beads
from ipi.utils.depend import *
from ipi.utils.units import *

def main(filename):

   ipos=open(filename,"r")
   
   natoms = 0
   ifr = 0
   while True:
      try:
         pos, cell = read_xyz(ipos, readcell=True)
         cell.array_pbc(pos.q)
      except EOFError: # finished reading files
         sys.exit(0)

      print_pdb(pos, cell)
      ifr+=1


if __name__ == '__main__':
   main(*sys.argv[1:])
