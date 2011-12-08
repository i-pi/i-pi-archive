import numpy as np
import math, sys
import utils.cell_convert
from engine.cell import *
from engine.atoms import *
from utils.units import *

def print_pdb(atoms, ncell, filedesc = sys.stdout):
   """Takes the system and gives pdb formatted output for the unit cell and the
      atomic positions """

   a, b, c, alpha, beta, gamma = utils.cell_convert.h2abc(ncell.h)
   alpha *= 180.0/math.pi #radian to degree conversion
   beta  *= 180.0/math.pi
   gamma *= 180.0/math.pi
   
   z = 1 #number of polymeric chains in a unit cell. I can't decide if 1 or 0 is more sensible for this...

   filedesc.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%s%4i\n" % (a, b, c, alpha, beta, gamma, " P 1        ", z))

   for i in range(0,len(atoms)): 
      filedesc.write("ATOM  %5i %4s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2i\n" % (i+1, atoms[i].name,' ','  1',' ',1,' ',atoms[i].q[0],atoms[i].q[1],atoms[i].q[2],0.0,0.0,'  ',0))

   filedesc.write("END\n")

def read_pdb(filedesc):
   """Takes a pdb-style file and creates a system with the appropriate unit
      cell and atom positions"""

   header = filedesc.readline()
   a = float(header[6:15]);      b = float(header[15:24]);    c = float(header[24:33]);
   alpha = float(header[33:40]); beta = float(header[40:47]); gamma = float(header[47:54]);
   alpha *= math.pi/180.0;       beta *= math.pi/180.0;       gamma *= math.pi/180.0
   h = utils.cell_convert.abc2h(a, b, c, alpha, beta, gamma)
   cell=Cell(h)
   
   natoms = 0
   body = filedesc.readline()
   qatoms=[]; names=[]
   while (body != "" and body != "END"):
      natoms += 1
      names.append(body[12:16].strip())
      x = float(body[31:39])
      y = float(body[39:47])
      z = float(body[47:55])
      pos = numpy.array([x,y,z])
      qatoms.append(pos)
      
      body = filedesc.readline()
   
   atoms=Atoms(natoms)
   for i in range(natoms):
      atoms[i].q=qatoms[i]
      atoms[i].name=names[i]
      atoms[i].m=Elements.mass(names[i])

   print atoms.q
   return atoms, cell
