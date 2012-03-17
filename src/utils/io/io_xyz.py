"""Contains the functions used to print the trajectories and read input
configurations.

Functions:
   print_xyz_path: Prints all the bead configurations, and shows the ring
      polymer connectivity.
   print_xyz: Prints the centroid configurations.
   read_xyz: Reads the cell parameters and atom configurations from a xyz file.
"""

__all__ = ['print_xyz_path', 'print_xyz', 'read_xyz']

import numpy as np
import math, sys
import utils.mathtools as mt
from engine.cell import Cell
from engine.atoms import Atoms
from utils.units import *
import gc

def print_xyz_path(beads, cell, filedesc = sys.stdout):
   """Prints all the bead configurations, into a xyz formatted file.

   Prints all the replicas for each time step separately, rather than all at
   once.

   Args:
      beads: A beads object giving the bead positions.
      cell: A cell object giving the system box.
      filedesc: An open writable file object. Defaults to standard output.
   """

   a, b, c, alpha, beta, gamma = mt.h2abc(cell.h)
   alpha *= 180.0/math.pi 
   beta  *= 180.0/math.pi
   gamma *= 180.0/math.pi

   natoms = beads.natoms
   nbeads = beads.nbeads   
   for j in range(nbeads):
      filedesc.write("%d\n# bead: %d CELL(abcABC): %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f \n" % ( natoms, j, a,b,c,alpha,beta,gamma))
      for i in range(natoms):
         bead = beads[j][i]
         filedesc.write("%8s %12.5e %12.5e %12.5e\n" % (bead.name[0],bead.q[0],bead.q[1],bead.q[2]))

def print_xyz(atoms, cell, filedesc = sys.stdout, title=""):
   """Prints the centroid configurations, into a xyz formatted file.

   Args:
      beads: An atoms object giving the centroid positions.
      cell: A cell object giving the system box.
      filedesc: An open writable file object. Defaults to standard output.
   """

   a, b, c, alpha, beta, gamma = mt.h2abc(cell.h)
   alpha *= 180.0/math.pi 
   beta  *= 180.0/math.pi
   gamma *= 180.0/math.pi

   natoms = atoms.natoms
   filedesc.write("%d\n# CELL(abcABC): %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %s\n" % ( natoms, a,b,c,alpha,beta,gamma, title))
   for i in range(natoms):
      atom = atoms[i]
      filedesc.write("%8s %12.5e %12.5e %12.5e\n" % (atom.name[0],atom.q[0],atom.q[1],atom.q[2]))   

def read_xyz(filedesc):
   """Takes a xyz-style file and creates an Atoms and Cell object.

   Args:
      filedesc: An open readable file object from a xyz formatted file.

   Returns:
      An Atoms object with the appropriate atom labels, masses and positions, 
      and a Cell object with the appropriate cell dimensions and an estimate 
      of a reasonable cell mass.
   """

   print "WILL DO THIS SOMETIMES...."

   return atoms, cell
