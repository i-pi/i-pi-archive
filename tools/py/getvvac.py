#!/usr/bin/env python2

""" getvvac.py

Computes velocity autocorrelation functions from i-pi outputs.

Assumes the input files are in xyz format and atomic units,
with prefix.vc.xyz naming scheme.

Syntax:
   getvvac.py prefix lag
"""


import numpy as np

import sys
from ipi.utils.io import iter_file_name, read_file
from ipi.utils.depend import *
from ipi.utils.units import *


def main(prefix, mlag, pad, label=None):
   # TODO: This really needs an argument parser to allow for several labels 
   mlag = int(mlag)
   npad = int(pad)

   ff = open(prefix + ".vc.xyz")
   rr = read_file("xyz", ff, output = "array")
   threenatoms = len(rr['data'])
   natoms = threenatoms / 3
   labelbool = np.ones(natoms, bool)
   if(label != None):
      labelbool = rr['names'] == label
   ff.close()

   ifile = open(prefix+".vc.xyz")
   ofile1 = prefix + "_" + "raw" + "_" + str(mlag) + ".vvft"
   ofile2 = prefix + "_" + "win" + "_" + str(mlag) + "_" + str(npad) + ".vvft"

   vel = np.zeros((2 * mlag, labelbool.sum(), 3) , float)
   fvvac = np.zeros((2 * mlag) / 2 + 1, float)
   omega = np.asarray(range(len(fvvac)))/float(len(fvvac))
   dt = 1.0 / float(2 * mlag)
   count = 0

   while True:

     try :
        for i in range(2*mlag):
            rr = read_file("xyz", ifile, output="array")
            vel[i] = rr['data'].reshape(natoms,3)[labelbool]

        fvel = np.fft.rfft(vel , axis = 0)
        tfvvac = fvel * np.conjugate(fvel)

        fvvac = fvvac + 3.0 * np.real(np.mean(tfvvac, axis = (1,2))) * dt / (2 * np.pi) # / (2 * mlag + npad) / 2 * mlag
        count = count + 1

     except EOFError:
        break

   fvvac = np.real(fvvac) / count 
   np.savetxt(ofile1, np.vstack((omega, fvvac)).T)

   vvac = np.fft.ihfft(fvvac) * np.bartlett(2 * mlag)[mlag:]
   pfvvac = np.real(np.fft.hfft(vvac , axis = 0, n = 2 * mlag + npad))[:2 * mlag]
   omega = (np.asarray(range(2 * mlag + npad))/float(2 * mlag + npad))[:2 * mlag]
   np.savetxt(ofile2, np.vstack((omega, pfvvac)).T[:mlag])
   

if __name__ == '__main__':
   main(*sys.argv[1:])
