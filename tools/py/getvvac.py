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

   ifile=open(prefix+".vc.xyz")
   ofile=prefix + "_" + str(mlag) + ".vvft"


   vel = np.zeros((2*mlag, labelbool.sum(), 3) , float)
   rvvac = np.zeros(2*mlag + npad, float)
   omega = np.asarray(range(2*mlag + npad) )/float(2*mlag + npad)
   time = np.asarray(range(mlag))
   window=np.bartlett(2*mlag).reshape((2*mlag,1,1))
   count = 0

   while True:

     try :
        for i in range(2*mlag):
            rr = read_file("xyz", ifile, output="array")
            vel[i] = rr['data'].reshape(natoms,3)[labelbool]

        vel = vel - np.mean(vel, axis=0)

        tmp = np.fft.fft(vel*window, n=2*mlag+npad , axis=0, norm="ortho")
        tmp = tmp * np.conjugate(tmp)
        tmp = np.real(np.mean(tmp, axis=(1,2)))
        rvvac = rvvac + tmp
        count = count + 1
        print count

     except EOFError:
        break
   
   np.savetxt(ofile, np.vstack((np.real(omega)[0:2*mlag], np.real(rvvac/count)[0:2*mlag] * (mlag + npad) / mlag )).T[0:mlag])

if __name__ == '__main__':
   main(*sys.argv[1:])
