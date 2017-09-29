#!/usr/bin/env python2 

description = """
   Computes the 3D momentum distribution having as input the end-to-end distances in atomic units,
   stored as a file with one sample per row, e.g. 
      Delta_x  Delta_y  Delta_z
      ....
   The output contains both the histogram of deltas 
   Dx Dy Dz n(D) 
   and the 3D PMD in atomic units, printed as a grid file
   px py pz n(p)   
"""

import argparse
import numpy as np

def histo3d(data, xg, yg, zg, sigma):
    ly=xg*0.0
    for x in data:
        ly += np.exp(-((xg-x[0])**2 + (yg-x[1])**2+ (zg-x[2])**2)*(0.5*sigma**2))
    return ly

def get_np(prefix, fname, bsize, P, m, Tkelv, nskip, s, ns):   

    T= 3.1668105e-06*Tkelv
    print "Kernel width: %10.5e - Grid spacing: %10.5e" % (1/np.sqrt(T * P * m), 2*s/(ns-1))
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
    pxgrid= pxgrid- pxstep/2.

    pyi = -np.pi/(dqygrid[1]-dqygrid[0])
    pyf = +np.pi/(dqygrid[1]-dqygrid[0])
    pystep = 2* np.pi / np.abs(dqygrid[-1]-dqygrid[0])
    pygrid = np.linspace(pyi,pyf,ns)
    pygrid= pygrid- pystep/2.

    pzi = -np.pi/(dqzgrid[1]-dqzgrid[0])
    pzf = +np.pi/(dqzgrid[1]-dqzgrid[0])
    pzstep = 2* np.pi / np.abs(dqzgrid[-1]-dqzgrid[0])
    pzgrid = np.linspace(pzi,pzf,ns)
    pzgrid= pzgrid- pzstep/2.
 
    pgrid= []
    for x in range(len(pxgrid)):
       for y in range(len(pygrid)):
          for z in range(len(pzgrid)):                  
             pgrid.append([pxgrid[x],pygrid[y],pzgrid[z]])

    pgrid= np.asarray(pgrid)


    #Read the end to end distances from file
    data_path =str(fname)
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
    
    np.savetxt(str(prefix + ".np"), np.c_[pgrid,avgnp3d])
 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-T", type=float, required=True, help="Specify the temperature of the system in Kelvin [required]")
    parser.add_argument("-P", type=int, required=True, help="Specify the number of beads [required]")
    parser.add_argument("-m", type=float, default= 1837 , help="Specify the mass of the atom in atomic units-default is hydrogen")
    parser.add_argument("-block", type=int, default=-4, help="Specify the size of the blocks [defaults: 1/4th of the trajectory]")
    parser.add_argument("--prefix",type=str, default="pmd3d", help="Prefix for output files")
    parser.add_argument("--fname",type=str,default="", help="name of the end-to-end distances file")    
    parser.add_argument("--skip", type=int, default= 50 , help="Removes the equilibration steps")
    parser.add_argument("-dint", type=float, default=15.0, help="Specify the positive extrema of the interval to build the histogram ([-dint,dint])")
    parser.add_argument("-ns", type=int, default=60, help="Specify the number of point to use for the histogram")
    args = parser.parse_args()

    get_np(args.prefix, args.fname, args.block, args.P, args.m, args.T, args.skip, args.dint, args.ns)


