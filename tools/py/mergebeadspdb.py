#!/usr/bin/env python2

""" mergebeadspdb.py

Reads positions of individual beads from an i-PI run and
assemles them in a pdb describing the ring polymer connectivity.

Assumes the input files are in pdb format names prefix.pos_*.pdb.

Syntax:
   mergebeadspdb.py prefix
"""


import numpy as np
import sys, glob
from ipi.utils.io import read_file, print_file_path
from ipi.engine.beads import Beads
from ipi.engine.cell import Cell
from ipi.utils.depend import *
from ipi.utils.units import *


def main(prefix):

   ipos=[]
   imode=[]
   for filename in sorted(glob.glob(prefix+".pos*")):
      imode.append(filename.split(".")[-1])
      ipos.append(open(filename,"r"))

   nbeads = len(ipos)
   natoms = 0
   ifr = 0
   while True:
      try:
         for i in range(nbeads):
            ret = read_file(imode[i], ipos[i], readcell="true")
            pos = ret["atoms"]
            cell = ret["cell"]
            if natoms == 0:
               natoms = pos.natoms
               beads = Beads(natoms,nbeads)
            beads[i].q = pos.q
            beads.names = pos.names
      except EOFError: # finished reading files
         sys.exit(0)

      print_file_path("pdb", beads, cell)
      ifr+=1


if __name__ == '__main__':
   main(*sys.argv[1:])
