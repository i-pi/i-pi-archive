#!/usr/bin/env python2 
"""

Relies on the infrastructure of i-pi, so the ipi package should
be installed in the Python module directory, or the i-pi
main directory must be added to the PYTHONPATH environment variable.

Cuts short the output of a previous i-pi simulation, up to the
step indicated in the <step> field of the input file.
This is useful to restart a simulation that crashed.

It should be run in the same dyrectory as where i-pi was (or is being)
run, and simply fetches all information from the simulation input file.
One should also specify a directory name in which the trimmed files
will be output.

Syntax:
   trimsim.py inputfile.xml
"""


import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from ipi.utils.io import read_file
from ipi.engine.outputs import *
from ipi.engine.properties import getkey
from ipi.inputs.simulation import InputSimulation
from ipi.utils.io.inputs import io_xml
from ipi.utils.units import unit_to_internal, unit_to_user


def kernel(x, mean=0, sigma=1):
    return np.exp(-(x- mean)**2*(0.5*sigma**2))

def histo(data, delta, k, mean, sigma):
    ly=delta*0.0
    for x in data:
        ly+=k(delta-x, mean, sigma)
    return ly

def get_np(path, fname, bsize, P, m, T, nskip, s, ns):   
    # initialises the data files.
    data_1 = np.zeros((bsize, 3) , float)
    data_2 = np.zeros((bsize, 3) , float)
    dq = np.zeros((bsize,3) , float)
    dqxgrid = np.linspace(-s, s, ns)
    dqygrid = np.linspace(-s, s, ns)
    dqzgrid = np.linspace(-s, s, ns)
    nplistx = []
    nplisty = []
    nplistz = []

    #Read the end to end distances from file
    data_path =str(path + fname)
    delta= np.loadtxt(data_path)
    step = np.shape(delta)[0] 
   
    n_block =int(step/bsize)
    if (n_block ==0):
             print 'not enough data to build a block'
             exit()
    for x in xrange(n_block):
        dq = delta[x*bsize : (x+1)*bsize]
        hx = histo(np.concatenate((dq.T[0], -dq.T[0])), dqxgrid, kernel, 0, np.sqrt(T * P * m))
        hy = histo(np.concatenate((dq.T[1], -dq.T[1])), dqygrid, kernel, 0, np.sqrt(T * P * m))
        hz = histo(np.concatenate((dq.T[2], -dq.T[2])), dqzgrid, kernel, 0, np.sqrt(T * P * m))
       
        # Defines the grid for momentum.
        pxi = -np.pi/(dqxgrid[1]-dqxgrid[0])
        pxf = +np.pi/(dqxgrid[1]-dqxgrid[0])
        pxstep = 2* np.pi / np.abs(dqxgrid[-1]-dqxgrid[0])
        pxgrid = np.linspace(pxi,pxf,ns)

	pyi = -np.pi/(dqygrid[1]-dqygrid[0])
        pyf = +np.pi/(dqygrid[1]-dqygrid[0])
        pystep = 2* np.pi / np.abs(dqygrid[-1]-dqygrid[0])
        pygrid = np.linspace(pyi,pyf,ns)

        pzi = -np.pi/(dqzgrid[1]-dqzgrid[0])
        pzf = +np.pi/(dqzgrid[1]-dqzgrid[0])
        pzstep = 2* np.pi / np.abs(dqzgrid[-1]-dqzgrid[0])
        pzgrid = np.linspace(pzi,pzf,ns)
            

        # Computes the Fourier transform of the end to end vector.
        npx = np.abs(np.fft.fftshift(np.fft.fft(hx)))
        npy = np.abs(np.fft.fftshift(np.fft.fft(hy)))
        npz = np.abs(np.fft.fftshift(np.fft.fft(hz)))
           
        nplistx.append(npx)
        nplisty.append(npy)
        nplistz.append(npz)
    
    avgnpx = np.mean(np.asarray(nplistx), axis = 0)
    avgnpy = np.mean(np.asarray(nplisty), axis = 0)
    avgnpz = np.mean(np.asarray(nplistz), axis = 0)
    normx=np.sum(avgnpx)
    normy=np.sum(avgnpy)
    normz=np.sum(avgnpz)
    errnpx = np.std(np.asarray(nplistx), axis = 0)/ np.sqrt(n_block)/normx
    avgnpx= avgnpx/normx
    errnpy = np.std(np.asarray(nplisty), axis = 0)/ np.sqrt(n_block)/normy
    avgnpy= avgnpy/normy
    errnpz = np.std(np.asarray(nplistz), axis = 0)/ np.sqrt(n_block)/normz
    avgnpz= avgnpz/normz

    avgpsqnpx = pxgrid**2*avgnpx/pxstep
    errpsqnpx = pxgrid**2*errnpx/pxstep
    avgpsqnpy = pygrid**2*avgnpy/pystep
    errpsqnpy = pygrid**2*errnpy/pystep
    avgpsqnpz = pzgrid**2*avgnpz/pzstep
    errpsqnpz = pzgrid**2*errnpz/pzstep
    
    np.savetxt(str(path + "np.data"), np.c_[pxgrid,avgnpx,errnpx,avgnpy,errnpy,avgnpz,errnpz])
    np.savetxt(str(path + "psq-np.data"), np.c_[pxgrid,avgpsqnpx,errpsqnpx,avgpsqnpy,errpsqnpy,avgpsqnpz,errpsqnpz])
    
    
    psqmedx =  0.
    psqmed2x = 0.
    psqmedy =  0.
    psqmed2y = 0.
    psqmedz =  0.
    psqmed2z = 0.
    for i in range(n_block):
         psqmedx= psqmedx + np.dot(pxgrid**2,np.asarray(nplistx)[i,:])/normx
         psqmed2x = psqmed2x + (np.dot(pxgrid**2,np.asarray(nplistx)[i,:])/normx)**2
         psqmedy= psqmedy + np.dot(pygrid**2,np.asarray(nplisty)[i,:])/normy
         psqmed2y = psqmed2y + (np.dot(pygrid**2,np.asarray(nplisty)[i,:])/normy)**2
         psqmedz= psqmedz + np.dot(pzgrid**2,np.asarray(nplistz)[i,:])/normz
         psqmed2z = psqmed2z + (np.dot(pzgrid**2,np.asarray(nplistz)[i,:])/normz)**2
         
    print 'number of blocks', n_block
    print 'av_px^2', psqmedx/n_block, 'sigmax', np.sqrt((psqmed2x/n_block) - (psqmedx/n_block)**2)/np.sqrt(n_block)
    print 'av_py^2', psqmedy/n_block, 'sigmay', np.sqrt((psqmed2y/n_block) - (psqmedy/n_block)**2)/np.sqrt(n_block)
    print 'av_pz^2', psqmedz/n_block, 'sigmaz', np.sqrt((psqmed2z/n_block) - (psqmedz/n_block)**2)/np.sqrt(n_block)
 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument("projpath")
    parser.add_argument("--fname",type=str,default="", help="name of the end-to-end distances file")
    parser.add_argument("-bsize", type=int, default=30000, help="Specify the size of the blocks")
    parser.add_argument("-P", type=int, default= 1 , help="Specify the number of beads")
    parser.add_argument("-m", type=float, default= 1837 , help="Specify the mass of the atom in atomic units")
    parser.add_argument("-T", type=float, default= 0.00095004315 , help="Specify the temperature of the system in hartree (kb*T)")
    parser.add_argument("-nskip", type=int, default= 100 , help="Removes the equilibration steps")
    parser.add_argument("-dint", type=float, default=15.0,help="Specify the positive extrema of the interval to build the histogram ([-dint,dint])")
    parser.add_argument("-ns", type=float, default=10000, help="Specify the number of point to use for the histogram")
    args = parser.parse_args()

    get_np(args.projpath, args.fname, args.bsize, args.P, args.m, args.T, args.nskip, args.dint, args.ns)


