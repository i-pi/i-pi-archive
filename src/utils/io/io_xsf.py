"""Contains the functions used to print the trajectories and read input
configurations with xsf formatting.

Functions:
   print_xsf: Prints the centroid configurations.
   read_xsf: Reads the cell parameters and atom configurations from a xsf file.
"""

#__all__ = ['print_xsf_path', 'print_xsf', 'read_xsf']
__all__ = ['print_xsf', 'read_xsf']

import numpy as np
import math, sys
from utils.depend import depstrip
import utils.mathtools as mt
from engine.atoms import Atoms
from engine.cell import Cell
from utils.units import *
import gc

#def print_xsf_path(beads, cell, filedesc = sys.stdout):
#   """Prints all the bead configurations, into a xsf formatted file.
#   Prints all the replicas for each time step separately, rather than all at
#   once.
#
#   Args:
#      beads: A beads object giving the bead positions.
#      cell: A cell object giving the system box.
#      filedesc: An open writable file object. Defaults to standard output.
#   """
#
#   a, b, c, alpha, beta, gamma = mt.h2abc(cell.h)
#   alpha *= 180.0/math.pi 
#   beta  *= 180.0/math.pi
#   gamma *= 180.0/math.pi
#
#   natoms = beads.natoms
#   nbeads = beads.nbeads   
#   for j in range(nbeads):
#      filedesc.write("%d\n# bead: %d CELL(abcABC): %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f \n" % ( natoms, j, a,b,c,alpha,beta,gamma))
#      for i in range(natoms):
#         bead = beads[j][i]
#         filedesc.write("%8s %12.5e %12.5e %12.5e\n" % (bead.name[0],bead.q[0],bead.q[1],bead.q[2]))

def print_xsf(atoms, cell, filedesc = sys.stdout, step = 0, total_steps = 1):
   """Prints the centroid configurations, into a xsf formatted file.

   Args:
      beads: An atoms object giving the centroid positions.
      cell: A cell object giving the system box.
      filedesc: An open writable file object. Defaults to standard output.
      step: The current step.
      total_steps: The total number of simulation steps.
   """

   natoms = atoms.natoms
   h = "%10.5f %10.5f %10.5f\n %10.5f %10.5f %10.5f\n %10.5f %10.5f %10.5f\n" % (cell.h[0,0], cell.h[1,0], cell.h[2,0], cell.h[0,1], cell.h[1,1], cell.h[2,1], cell.h[0,2], cell.h[1,2], cell.h[2,2])

   if step == 0:
      filedesc.write("ANIMSTEPS " + str(total_steps))
      filedesc.write("\nCRYSTAL\n")

   filedesc.write("PRIMVEC " + str(step+1))
   filedesc.write("\n" + h)
   filedesc.write("CONVVEC " + str(step+1))
   filedesc.write("\n" + h)
   filedesc.write("PRIMCOORD " + str(step+1) + "\n" + str(natoms) + " 1\n")

   # direct access to avoid unnecessary slow-down
   qs=depstrip(atoms.q)
   lab=depstrip(atoms.names)
   for i in range(natoms):
      filedesc.write("%8s %12.5e %12.5e %12.5e\n" % (lab[i],qs[3*i],qs[3*i+1],qs[3*i+2]))   
   filedesc.flush()

def read_xsf(filedesc):
   """Takes a xsf-style file and creates an Atoms object.

   Args:
      filedesc: An open readable file object from a xsf formatted file.

   Returns:
      An Atoms object with the appropriate atom labels, masses and positions.
   """

   ignore = False

   line = filedesc.readline().split()
   if line[0] == "ANIMSTEPS":
      print ".axsf file given as input, only the first frame will be used."
      line = filedesc.readline().split()

   if line[0] == "MOLECULE" or line[0] == "POLYMER" or line[0] == "SLAB":
      print "No cell specified, please specify it in xml input file."
      ignore = True
      line = filedesc.readline().split()
   elif line[0] != "CRYSTAL":
      print "No cell specified, please specify it in xml input file."
      ignore = True
   else:
      line = filedesc.readline().split() 

   qatoms = []
   names = []
   if not ignore:
      h = np.zeros((3,3))

   while (True):
      print line
      if (ignore and (line[0] != "PRIMCOORD" and line[0] != "ATOMS")):
         line = filedesc.readline().split()
      elif (line[0] == "PRIMVEC"):
         line = filedesc.readline().split()
         for i in range(3):
            x = float(line[0])
            y = float(line[1])
            z = float(line[2])
            h[:,i] = np.array([x,y,z])
            line = filedesc.readline().split()
      elif (line[0] == "CONVVEC"):
         line = filedesc.readline().split()
         ch = np.zeros((3,3))
         for i in range(3):
            x = float(line[0])
            y = float(line[1])
            z = float(line[2])
            ch[:,i] = np.array([x,y,z])
            line = filedesc.readline().split()
         if (h != ch).all():
            print "PRIMVEC != CONVVEC not implemented"
            exit()
      elif (line[0] == "PRIMCOORD"):
         line = filedesc.readline().split()
         natoms = int(line[0])
         line = filedesc.readline().split()
         while (len(line) == 4):
            names.append(line[0])
            x = float(line[1])
            y = float(line[2])
            z = float(line[3])
            pos = np.array([x,y,z])
            qatoms.append(pos)
            line = filedesc.readline().split()
         if (len(names) != natoms):
            print "Wrong number of atoms given in PRIMCOORD"
            exit()
         break
      elif (line[0] == "ATOMS"):
         line = filedesc.readline().split()
         while (len(line) == 4):
            names.append(line[0])
            x = float(line[1])
            y = float(line[2])
            z = float(line[3])
            pos = np.array([x,y,z])
            qatoms.append(pos)
            line = filedesc.readline().split()
         natoms = len(names)
         break
      elif (line[0] == "CONVCOORD"):
         print "CONVCOORD not currently implemented"
         exit()
      else:
         print "Unrecognized keyword, " + line[0] + ", in input file."
         exit()

   atoms = Atoms(natoms)
   totmass = 0.0
   for i in range(natoms):
      nat = atoms[i]
      nat.q = qatoms[i]
      nat.name = names[i]
      nat.m = Elements.mass(names[i])
      totmass += Elements.mass(names[i])

   if not ignore:
      cell = Cell(h)
      cell.m = totmass
   else:
      cell = None

   return atoms, cell
