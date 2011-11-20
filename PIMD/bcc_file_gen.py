from numpy import *
from engine import *
import sys
from engine import io_system
from engine import dynamics
from engine import forces
from engine import rp_engine
from engine import rp_dynamics
import numpy, random
from utils.mathtools import *

f = open("./testfile5.txt","r")
#must have 2*cell_no**3 atoms in it, with the appropriate names and cell

cell_no = 4
cell_len = 11.2423
syst = engine.System.from_pdbfile(f)
count = 0
for i in range(cell_no):
   for j in range(cell_no):
      for k in range(cell_no):
         syst.atoms[count].q.get_array()[:] = numpy.array([i*cell_len, j*cell_len, k*cell_len])
         count += 1
         syst.atoms[count].q.get_array()[:] = numpy.array([(i+0.5)*cell_len, (j+0.5)*cell_len, (k+0.5)*cell_len])
         count += 1
   
g = open("./bcc%1icell.pdb" % (cell_no), "w")
io_system.print_pdb(syst.atoms, syst.cell, g)
