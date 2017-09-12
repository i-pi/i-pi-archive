#!/usr/bin/env python2 
"""

Computes the momentum distribution having as input the end-to-end distances in atomic units.
The results are in atomic units.

"""

import argparse
import numpy as np


def histo3d(data,xg, yg, zg, sigma):
    ly=xg*0.0
    for x in data:
        ly+= np.exp(-((xg-x[0])**2 + (yg-x[1])**2+ (zg-x[2])**2)*(0.5*sigma**2))
    return ly

def get_np(path, fname, bsize, P, m, Tkelv, nskip, s, ns):   

    T= 3.1668105e-06*Tkelv
    dq = np.zeros((bsize,3) , float)
    dqxgrid = np.linspace(-s, s, ns)
    dqygrid = np.linspace(-s, s, ns)
    dqzgrid = np.linspace(-s, s, ns)
    #make a meshgrid for points
    xxx, yyy, zzz = np.meshgrid(dqxgrid, dqygrid, dqzgrid)


    # Defines the final grid for momentum.
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
 
    pgrid= []
    for x in range(len(pxgrid)):
       for y in range(len(pygrid)):
          for z in range(len(pzgrid)):                  
             pgrid.append([pxgrid[x],pygrid[y],pzgrid[z]])

    pgrid= np.asarray(pgrid)


    #Read the end to end distances from file
    data_path =str(path + fname)
    delta= np.loadtxt(data_path)
    step = np.shape(delta)[0] 
    
    nplist3d=[]
    n_block =int(step/bsize)
    if (n_block ==0):
             print 'not enough data to build a block'
             exit()
    for x in xrange(n_block):
        dq = delta[x*bsize : (x+1)*bsize]
        h3d = histo3d(np.concatenate((dq, -dq)), xxx, yyy, zzz, np.sqrt(T * P * m)) 
               
        # Computes the 3D Fourier transform of the histogram of the end-to-end distances.
        npvec3d= np.abs(np.fft.fftshift(np.fft.fftn(h3d)))       

        nplist3d.append(np.reshape(npvec3d,(ns*ns*ns)))
        
    avgnp3d= np.mean(np.asarray(nplist3d), axis = 0)
    norm=np.sum(avgnp3d)
    avgnp3d= avgnp3d/norm
    
    np.savetxt(str(path + "np3d.data"), np.c_[pgrid,avgnp3d])
 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument("--path",type=str, default="", help="path of the folder conatining the end-to-end distances file")
    parser.add_argument("--fname",type=str,default="", help="name of the end-to-end distances file")
    parser.add_argument("-bsize", type=int, default=10000, help="Specify the size of the blocks")
    parser.add_argument("-P", type=int, default= 1 , help="Specify the number of beads")
    parser.add_argument("-m", type=float, default= 1837 , help="Specify the mass of the atom in atomic units-default is hydorgen")
    parser.add_argument("-T", type=float, default= 300 , help="Specify the temperature of the system in kelvin")
    parser.add_argument("-nskip", type=int, default= 50 , help="Removes the equilibration steps")
    parser.add_argument("-dint", type=float, default=15.0, help="Specify the positive extrema of the interval to build the histogram ([-dint,dint])")
    parser.add_argument("-ns", type=int, default=60, help="Specify the number of point to use for the histogram")
    args = parser.parse_args()

    get_np(args.path, args.fname, args.bsize, args.P, args.m, args.T, args.nskip, args.dint, args.ns)


