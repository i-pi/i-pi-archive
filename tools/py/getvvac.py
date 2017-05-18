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

   ff = open(prefix+".vc.xyz")
   rr = read_file("xyz", ff, output="array")
   threenatoms = len(rr['data'])
   natoms = threenatoms/3
   labelbool = np.ones(natoms, bool)
   if(label != None):
      labelbool = rr['names'] == label
   ff.close()

   ifile = open(prefix+".vc.xyz")
   ofile1 = prefix + "_" + "raw" + "_" + str(mlag) + ".vvft"
   ofile2 = prefix + "_" + "win" + "_" + str(mlag) + "_" + str(npad) + ".vvft"


   vel = np.zeros((2*mlag, labelbool.sum(), 3) , float)
   fvvac = np.zeros(2*mlag, float)
   omega = np.asarray(range(2*mlag))/float(2*mlag)
   time = np.asarray(range(mlag))
   window = np.bartlett(2*mlag)[mlag:]
   count = 0

   while True:

     try :
        for i in range(2*mlag):
            rr = read_file("xyz", ifile, output="array")
            vel[i] = rr['data'].reshape(natoms,3)[labelbool]

        vel = vel - np.mean(vel, axis=0)

        tmp = np.fft.fft(vel, axis=0, norm="ortho") * 2.0/np.sqrt(2.0*np.pi)
        tmp = tmp * np.conjugate(tmp)

        fvvac = fvvac + 3.0 * np.real(np.mean(tmp, axis=(1,2)))
        count = count + 1

     except EOFError:
        break

   fvvac = np.real(fvvac / count * np.sqrt(2.0*mlag))
   np.savetxt(ofile1, np.vstack((omega, fvvac)).T )

   vvac = np.fft.ihfft(fvvac, norm="ortho")[:mlag]
   wvvac = np.array([i*j for i,j in zip(window, vvac)])
   fwvvac = np.fft.hfft(wvvac, n=2*mlag+npad, norm="ortho")  / np.sqrt(2.0*mlag) * np.sqrt(2.0*mlag + npad)
   omega = np.asarray(range(2*mlag+npad))/float(2*mlag+npad)

   np.savetxt(ofile2, np.vstack((omega, np.real(fwvvac))).T[0:2*mlag])

if __name__ == '__main__':
   main(*sys.argv[1:])
