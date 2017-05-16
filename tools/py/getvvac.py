#!*/usr/bin/env python2

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


def main(prefix, mlag, label):
   # TODO: This really needs an argument parser to allow for several labels 
   mlag = int(mlag)


   ff = open(prefix+".vc.xyz")
   rr = read_file("xyz", ff, output="array")
   threenatoms = len(rr['data'])
   natoms = threenatoms/3
   labelbool = rr['names'] == label
   # testhack
   #a=rr['names'] == "O"
   #b=rr['names'] == "H"
   #labelbool=a+b
   # end testhack
   ff.close()

   ifile=open(prefix+".vc.xyz")
   ofile=prefix+"_"+str(mlag)+".vvft"
   ofile2=prefix+"_"+str(mlag)+"_window.vvft"
   ofile3=prefix+"_"+str(mlag)+".vvac"


   vel = np.zeros((2*mlag, labelbool.sum(), 3) , float)
   rvvac = np.zeros(2*mlag, float)
   omega = np.asarray(range(2*mlag))/float(2*mlag)
   time = np.asarray(range(mlag))
   count = 0
# TODO: pad should become an argument as well
   pad=10000

   while True:

     try :
        for i in range(2*mlag):
            rr = read_file("xyz", ifile, output="array")
            vel[i] = rr['data'].reshape(natoms,3)[labelbool]

        vel = vel - np.mean(vel, axis=0)

        tmp = np.fft.fft(vel, axis=0, norm="ortho")
        tmp = tmp * np.conjugate(tmp)
        tmp = np.real(np.mean(tmp, axis=(1,2)))
        rvvac = rvvac + tmp
        count = count + 1

     except EOFError:
        break

# Go for windowing and padding
# First the autocorrelation function. "ortho" gives back the 1/sqrt(mlag) factor
   autocorr=np.fft.ihfft(rvvac, norm="ortho")[:mlag]
# Let's try first the triangular window
   window=np.bartlett(2.*mlag)[mlag:]
# Now multiply them in real time (could we do it in Fourier space?)
   autowindow=np.array([i*j for i,j in zip(window, autocorr)])
# And finally, re-fourier transform assuming the signal is hermitian in order to avoid noise
# Giving n greater than the array dimension pads it with zeroes
   windowvvac=np.fft.hfft(autowindow, n=2*mlag+pad, norm="ortho")

   omegapad = np.asarray(range(2*mlag+pad))/float(2*mlag+pad)

   np.savetxt(ofile3, np.vstack((np.real(time), np.real(autocorr), np.real(autowindow[:mlag]))).T[0:2*mlag])
 
   np.savetxt(ofile, np.vstack((np.real(omega), np.real(rvvac/count))).T[0:mlag])

# normalize to maintain same area as before padding and account for the various 2*pi factors
   np.savetxt(ofile2, np.vstack((np.real(omegapad), np.real(windowvvac*(mlag+pad)/(count*2.*np.pi*mlag)))).T[0:2*mlag])

if __name__ == '__main__':
   main(*sys.argv[1:])
