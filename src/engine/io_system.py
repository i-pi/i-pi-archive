import numpy
import math
import cell 

def print_pdb(atoms, ncell):
   a, b, c, alpha, beta, gamma = cell.h2abc(ncell.h)
   alpha *= 180.0/math.pi
   beta  *= 180.0/math.pi
   gamma *= 180.0/math.pi
   
   z = 1 #number of polymeric chains in a unit cell. I can't decide if 1 or 0 is more sensible for this...
   print "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%s%4i" % (a, b, c, alpha, beta, gamma, " P 1        ", z)
   for i in range(0,len(atoms)): 
      print "ATOM  %5i %4s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2i" % (i+1, atoms[i].name,' ','  1',' ',1,' ',atoms[i].q[0],atoms[i].q[1],atoms[i].q[2],0.0,0.0,'  ',0)

def read_pdb(filedesc):

#We need a way of using these values to initialise the system

   header = filedesc.readline()
   a = float(header[6:15]);      b = float(header[15:24]);    c = float(header[24:33]);
   alpha = float(header[33:40]); beta = float(header[40:47]); gamma = float(header[47:54]);
   alpha *= math.pi/180.0;       beta *= math.pi/180.0;       gamma *= math.pi/180.0
   cell = numpy.array([a, b, c, alpha, beta, gamma])

   atoms = []
   natoms = 0
   body = filedesc.readline()
   while body != '':
      natoms += 1
      name = body[12:16]
      x = float(body[31:39])
      y = float(body[39:47])
      z = float(body[47:55])
      pos = numpy.array([x, y, z])
      atoms.append([name, pos]) 
      body = filedesc.readline()
   return atoms, cell, natoms
