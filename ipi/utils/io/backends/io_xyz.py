"""Functions used to read input configurations and print trajectories
in the XYZ format.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import math
import sys
import re

import numpy as np

import ipi.utils.mathtools as mt
from ipi.utils.depend import depstrip
from ipi.engine.atoms import Atoms
from ipi.engine.cell import Cell
from ipi.utils.units import *


__all__ = ['print_xyz_path', 'print_xyz', 'read_xyz', 'iter_xyz']


def print_xyz_path(beads, cell, filedesc = sys.stdout):
   """Prints all the bead configurations into a XYZ formatted file.

   Prints all the replicas for each time step separately, rather than all at
   once.

   Args:
      beads: A beads object giving the bead positions.
      cell: A cell object giving the system box.
      filedesc: An open writable file object. Defaults to standard output.
   """

   a, b, c, alpha, beta, gamma = mt.h2abc_deg(cell.h)

   natoms = beads.natoms
   nbeads = beads.nbeads
   for j in range(nbeads):
      filedesc.write("%d\n# bead: %d CELL(abcABC): %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f \n" % (natoms, j, a, b, c, alpha, beta, gamma))
      for i in range(natoms):
         qs = depstrip(beads.q)
         lab = depstrip(beads.names)
         filedesc.write("%8s %12.5e %12.5e %12.5e\n" % (lab[i], qs[j][3*i], qs[j][3*i+1], qs[j][3*i+2]))

def print_xyz(atoms, cell, filedesc = sys.stdout, title=""):
   """Prints an atomic configuration into an XYZ formatted file.

   Args:
      atoms: An atoms object giving the centroid positions.
      cell: A cell object giving the system box.
      filedesc: An open writable file object. Defaults to standard output.
      title: This gives a string to be appended to the comment line.
   """

   a, b, c, alpha, beta, gamma = mt.h2abc_deg(cell.h)

   natoms = atoms.natoms
   filedesc.write("%d\n# CELL(abcABC): %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %s\n" % ( natoms, a, b, c, alpha, beta, gamma, title))
   # direct access to avoid unnecessary slow-down
   qs = depstrip(atoms.q)
   lab = depstrip(atoms.names)
   for i in range(natoms):
      filedesc.write("%8s %12.5e %12.5e %12.5e\n" % (lab[i], qs[3*i], qs[3*i+1], qs[3*i+2]))

def read_xyz(filedesc, readcell=False):
   """Readss an XYZ-style file with i-pi style comments and creates an Atoms and Cell object

   Args:
      filedesc: An open readable file object from a xyz formatted file with i-PI header comments.

   Returns:
      An Atoms object with the appropriate atom labels, masses and positions.
      A Cell object.
   """

   natoms = filedesc.readline()
   if natoms == "":
      raise EOFError("The file descriptor hit EOF.")
   natoms = int(natoms)
   comment = filedesc.readline()
   recell = re.compile('# CELL.abcABC.: (.*) Traj')
   reres = recell.search(comment)
   if (reres is None):  # defaults to unit box
      h = mt.abc2h(1., 1., 1., np.pi/2, np.pi/2, np.pi/2)
   else:
      (a,b,c,alpha,beta,gamma) = reres.group(1).split()
      a = float(a)
      b = float(b)
      c = float(c)
      alpha =  float(alpha) * np.pi/180
      beta =  float(beta) * np.pi/180
      gamma =  float(gamma) * np.pi/180
      h = mt.abc2h(a, b, c, alpha, beta, gamma)
   cell = Cell(h)

   qatoms = []
   names = []
   masses = []
   iat = 0
   while (iat < natoms):
      body = filedesc.readline()
      if body.strip() == "":
         break
      body = body.split()
      name = body[0]
      names.append(name)
      masses.append(Elements.mass(name))
      x = float(body[1])
      y = float(body[2])
      z = float(body[3])
      qatoms.append(x)
      qatoms.append(y)
      qatoms.append(z)
      iat += 1

   if natoms != len(names):
      raise ValueError("The number of atom records does not match the header of the xyz file.")

   atoms = Atoms(natoms)
   atoms.q = np.asarray(qatoms)
   atoms.names = np.asarray(names, dtype='|S4')
   atoms.m = np.asarray(masses)

   if readcell: return atoms, cell
   else: return atoms

def iter_xyz(filedesc):
   """Takes a xyz-style file and yields one Atoms object after another.

   Args:
      filedesc: An open readable file object from a xyz formatted file.

   Returns:
      Generator over the xyz trajectory, that yields
      Atoms objects with the appropriate atom labels, masses and positions.
   """

   try:
      while 1:
         atoms = read_xyz(filedesc)
         yield atoms
   except EOFError:
      pass
