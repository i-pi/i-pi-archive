#!/usr/bin/env python2

""" getvvac.py
Computes velocity autocorrelation functions from i-pi outputs.
Assumes the input files are in xyz format and atomic units.
"""


import numpy as np
import argparse
import sys
from ipi.utils.io import iter_file_name, read_file
from ipi.utils.depend import *
from ipi.utils.units import *


def main():
   # adds description of the program.
    parser=argparse.ArgumentParser(description="Given the velocity of a system, computes the velocity-velocity autocorrelation function and its Fourier transform")

   # adds arguments.
    parser.add_argument("-ifile", "--input_file", required=True, nargs=1, type=str, default=None, help="the relative path to the xyz formatted velocity file")
    parser.add_argument("-mlag", "--maximum_lag", required=True, nargs=1, type=int, default=None, help="the maximum time lag for the autocorrelation function")
    parser.add_argument("-bsize", "--size_block", nargs=1, type=int, default=None,  help="the size of the block for ``chunk-by-chunk`` input.")
    parser.add_argument("-ftpad", "--size_padding", nargs=1, type=int, default=0, help="number of zeroes to be padded before the Fourier transform.")
    parser.add_argument("-ftbox", "--windowing", nargs=1, type=bool, default=True, help="if autocorrelation function should be zeroed at the boundaries before the Fourier transform.")
    parser.add_argument("-label", "--label", nargs=1, type=int, default=None, help="atomic species to be monitored")
    parser.add_argument("-oprefix", "--output_prefix", required=True, nargs=1, type=str, help="the prefix of the output file.")

    # parses the arguments
    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit()

    # stores the arguments
    ifile = str(args.input_file[-1])
    oprefix = str(args.output_prefix[-1])
    mlag = int(args.maximum_lag[-1])
    npad = int(args.size_padding[-1])
    ftbox = bool(args.size_padding[-1])
    if(args.size_block != None):
        bsize = int(args.size_block[-1])
        if (bsize < 2 * npad):
            raise ValueError("SIZE_BLOCK must be greater than 2 * MAXIMUM_LAG")
    else:
        bsize = npad
    if(args.label != None):
        label = str(args.label[-1])
    else:
        label = None

    #stores the number of atoms, and the index of the "chosen" atoms.
    ff = open(ifile)
    rr = read_file("xyz", ff, output = "array")
    threenatoms = len(rr['data'])
    natoms = threenatoms / 3
    labelbool = np.ones(natoms, bool)
    if(label != None):
       labelbool = rr['names'] == label
    ff.close()

    #initializes variables.
    ff = open(ifile)
    ofile = oprefix + ".vvft"
    vel = np.zeros((2 * bsize, labelbool.sum(), 3) , float)
    fvvac = np.zeros((2 * bsize) / 2 + 1, float)
    omega = (2 * np.pi) * np.asarray(range(2 * bsize)) / float(2 * bsize)
    dt = 1.0 / float(2 * bsize)
    count = 0
    if(ftbox == True):
        win = np.bartlett(2 * mlag + 1)
    else:
        win = np.ones(2 * mlag + 1, float)

    while True:

     try :
        #Reads the velocity in blocks.
        for i in range(2 * bsize):
            rr = read_file("xyz", ff, output="array")
            vel[i] = rr['data'].reshape(natoms,3)[labelbool]

        #Computes the Fourier transform of the velocity.
        fvel = np.fft.rfft(vel , axis = 0)

        #Computes the Fourier transform of the vvac applying the convolution theorem.
        tfvvac = fvel * np.conjugate(fvel)

        #Averages over all species and sums over the x,y,z directions. Also multiplies with the time step and a prefactor of (2pi)^-1.
        fvvac = fvvac + 3.0 * np.real(np.mean(tfvvac, axis = (1,2))) * dt / (2 * np.pi)
        count = count + 1

     except EOFError:
        break

    #Performs the block average of the Fourier transform.
    fvvac = np.real(fvvac) / count

    #Computes the inverse Fourier transform to get the vvac.
    vvac = np.fft.irfft(fvvac)[:mlag + 1]
    time = np.asarray(range(len(vvac))) * dt
    np.savetxt(oprefix + "-vvac.data" , np.vstack((time, vvac)).T[:mlag + npad])

    #Applies window in one direction and pads the vvac with zeroes.
    pvvac = np.append(vvac * win[mlag:], np.zeros(npad))

    #Recomputes the Fourier transform assuming the data is an even function of time.
    fpvvac = np.fft.hfft(pvvac)
    omega = (np.asarray(range(2 * (mlag + npad)))/float(2 * mlag + 2 * npad)) * (2 * np.pi)
    np.savetxt(oprefix + "-fvvac.data" , np.vstack((omega, fpvvac)).T[:mlag + npad])

if __name__ == '__main__':
   main()
