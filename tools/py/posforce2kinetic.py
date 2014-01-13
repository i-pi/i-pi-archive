#!/usr/bin/python
""" posforce2kinetic.py

Reads positions and forces from an i-PI run and computes the
centroid-virial kinetic energy estimator for each particle and
time frame, printing it out in the same format as it would have
been obtained had the run included the kinetic_cv and kinetic_od
output trajectories.

Assumes the input files are in xyz format and atomic units,
with prefix.pos_*.xyz and prefix.for_*.xyz naming scheme.

Syntax:
   trimsim.py prefix
"""

import numpy as np
import sys, glob
from ipi.utils.io.io_xyz import *
from ipi.engine.beads import Beads
from ipi.utils.depend import *

def main(prefix):

   ipos=[]
   for filename in glob.glob(prefix+".pos*"):
      ipos.append(open(filename,"r"))

   ifor=[]
   for filename in glob.glob(prefix+".for*"):
      ifor.append(open(filename,"r"))

   ikin=open(prefix+".kin.xyz","w")
   ikod=open(prefix+".kod.xyz","w")

   nbeads = len(ipos)
   if (nbeads!=len(ifor)): raise ValueError("Mismatch between number of output files for forces and positions")
   natoms = 0
   ifr = 1
   while True:
      try:
         for i in range(nbeads):
            pos = read_xyz(ipos[i])
            force = read_xyz(ifor[i])
            if natoms == 0:
               natoms = pos.natoms
               beads = Beads(natoms,nbeads)
               forces = Beads(natoms,nbeads)
               kcv = np.zeros((natoms,6),float)
            beads[i].q = pos.q
            forces[i].q = force.q
      except:
         raise
         sys.exit(0)

      q = depstrip(beads.q)
      f = depstrip(forces.q)
      qc = depstrip(beads.qc)
      kcv[:]=0
      for i in range(natoms):
         for j in range(nbeads):
            kcv[i,0] -= (q[j,i*3+0]-qc[i*3+0])*f[j,i*3+0]
            kcv[i,1] -= (q[j,i*3+1]-qc[i*3+1])*f[j,i*3+1]
            kcv[i,2] -= (q[j,i*3+2]-qc[i*3+2])*f[j,i*3+2]
            kcv[i,3] -= (q[j,i*3+0]-qc[i*3+0])*f[j,i*3+1]
            kcv[i,4] -= (q[j,i*3+0]-qc[i*3+0])*f[j,i*3+2]
            kcv[i,5] -= (q[j,i*3+1]-qc[i*3+1])*f[j,i*3+2]

      ifr+=1
      print "frame", ifr


if __name__ == '__main__':
   main(*sys.argv[1:])
