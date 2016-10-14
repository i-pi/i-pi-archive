#!/usr/bin/env python2

""" mergebeadspdb.py

Reads positions of individual beads from an i-PI run and
assemles them in a pdb describing the ring polymer connectivity.

Assumes the input files are in pdb format names prefix.pos_*.pdb.

Syntax:
   mergebeadspdb.py prefix
"""


import sys
from ipi.utils.io import read_file, print_file
from ipi.utils.depend import *
from ipi.utils.units import *


def main(filename):

   ipos=open(filename,"r")

   natoms = 0
   ifr = 0
   while True:
      try:
         ret = read_file("xyz", ipos, readcell=True)
         pos = ret["atoms"]
         cell = ret["cell"]
         cell.array_pbc(pos.q)
      except EOFError: # finished reading files
         sys.exit(0)

      print_file("pdb", pos, cell)
      ifr+=1


if __name__ == '__main__':
   main(*sys.argv[1:])
