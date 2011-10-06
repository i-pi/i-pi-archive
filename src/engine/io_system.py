import numpy
import math

def print_pdb(atoms, cell):
   a, b, c, alpha, beta, gamma = cell.h2abc()
   alpha *= 180.0/math.pi
   beta  *= 180.0/math.pi
   gamma *= 180.0/math.pi
   
   z = 1 #we need to find out what Z actually is
   print "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%s%4i" % (a, b, c, alpha, beta, gamma, " P 1       ", z)
   for i in range(0,len(atoms)): 
      print "ATOM  %5i %4s%1s%3s %1s%4i%1s%8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2i" % (i+1, atoms[i].name,'X','  1',' ',1,' ',atoms[i].q[0],atoms[i].q[1],atoms[i].q[2],0.0,0.0,'  ',0)

def read_pdb(filedesc):

#We need a way of using these values to initialise the system

   header = filedesc.readline()
   a = float(header[6:15]);      b = float(header[15:24]);     c = float(header[24:33]);
   alpha = float(header[33:40]); beta = float(header[40:47]);  gamma = float(header[47:54]);

   body = filedesc.readline()
   while body != '':
      name = body[12:16]
      x = float(body[31:39])
      y = float(body[39:47])
      z = float(body[47:55])
      pos = numpy.array([x, y, z])
      body = filedesc.readline()
   atoms = ''
   cell = ''
   return atoms, cell
