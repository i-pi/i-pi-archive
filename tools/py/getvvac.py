#!/usr/bin/env python2

""" kinetic2tag.py

Computes the Transient Anisotropic Gaussian (TAG) approximation
of the instantaneous kinetic energy tensor, with a moving average
triangular window of the specified lag. Needs files with
the per-atom diagonal and off-diagonal components of the kinetic
energy tensor estimator.

Assumes the input files are in xyz format and atomic units,
with prefix.kin.xyz and prefix.kod.xyz naming scheme.

Syntax:
   kinetic2tag.py prefix lag
"""


import numpy as np

import sys
from ipi.utils.io import iter_file_name, read_file
from ipi.utils.depend import *
from ipi.utils.units import *


def main(prefix, mlag, label):

   mlag = int(mlag)


   ff = open(prefix+".vc.xyz")
   rr = read_file("xyz", ff, output="array")
   threenatoms = len(rr['data'])
   natoms = threenatoms/3
   labelbool = rr['names'] == label
   ff.close()

   ifile=open(prefix+".vc.xyz")
   ofile=prefix+"_"+str(mlag)+".vvac"



   vel = np.zeros((2*mlag, labelbool.sum(), 3) , float)
   rvvac = np.zeros(2*mlag, float)
   omega = np.asarray(range(2*mlag))/float(2*mlag)
   count = 0

   while True:

     try :
        for i in range(2*mlag):
            rr = read_file("xyz", ifile, output="array")
            vel[i] = rr['data'].reshape(natoms,3)[labelbool]

        vel = vel - np.mean(vel, axis=0)

        tmp = np.fft.fft(vel, axis=0)
        tmp = tmp * np.conjugate(tmp)
        tmp = np.real(np.mean(tmp, axis=(1,2)))
        rvvac = rvvac + tmp
        count = count + 1

     except EOFError:
        break
 
   np.savetxt(ofile, np.vstack((omega, rvvac/count)).T[0:mlag])

if __name__ == '__main__':
   main(*sys.argv[1:])
