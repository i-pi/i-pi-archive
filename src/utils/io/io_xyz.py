"""Contains the functions used to print the trajectories and read input
configurations with xyz formatting.

Functions:
   print_xyz_path: Prints all the bead configurations.
   print_xyz: Prints the centroid configurations.
   read_xyz: Reads the cell parameters and atom configurations from a xyz file.
"""

__all__ = ['print_xyz_path', 'print_xyz', 'read_xyz']

import numpy as np
import math, sys
from utils.depend import depstrip
import utils.mathtools as mt
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
      title: This gives a string to be appended to the comment line.
   """

   a, b, c, alpha, beta, gamma = mt.h2abc(cell.h)
   alpha *= 180.0/math.pi
   beta  *= 180.0/math.pi
   gamma *= 180.0/math.pi

   natoms = atoms.natoms
   filedesc.write("%d\n# CELL(abcABC): %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %s\n" % ( natoms, a,b,c,alpha,beta,gamma, title))
   # direct access to avoid unnecessary slow-down
   qs=depstrip(atoms.q)
   lab=depstrip(atoms.names)
   for i in range(natoms):
      filedesc.write("%8s %12.5e %12.5e %12.5e\n" % (lab[i],qs[3*i],qs[3*i+1],qs[3*i+2]))

def read_xyz(filedesc):
   """Takes a xyz-style file and creates an Atoms object.

   Args:
      filedesc: An open readable file object from a xyz formatted file.

   Returns:
      An Atoms object with the appropriate atom labels, masses and positions.
   """

   natoms = int(filedesc.readline())
   comment = filedesc.readline()

   qatoms = []
   names = []
   iat=0
   while (iat<natoms):
      body = filedesc.readline()
      if body.strip() == "": break
      body = body.split()
      names.append(body[0])
      x = float(body[1])
      y = float(body[2])
      z = float(body[3])
      pos = np.array([x,y,z])
      qatoms.append(pos)
      iat+=1

   if natoms != len(names):
      raise ValueError("The number of atom records does not match the header of the xyz file.")

   atoms = Atoms(natoms)
   for i in range(natoms):
      nat = atoms[i]
      nat.q = qatoms[i]
      nat.name = names[i]
      nat.m = Elements.mass(names[i])

   return atoms
