#!/usr/bin/env python2 
description = """
Computes the quantum momentum distribution of a particle given the end-to-end distances.  
It computes both the three components of the momentum distribution and the spherically-averaged 
distribution of the proton momentum. It also computes <p^2> in the various directions, and the 
total contribute.
"""

import argparse
import numpy as np

def rad_kernel(x, delta, spread):
    # computes the kernel that describes the contribution to the radial average n(Delta) for a given spread factor 1/ispread
    if (x <= 1.0e-4): # use a Taylor expansion to get a stable result for small x
        delta2 = delta**2
        res= np.exp(-spread*delta2)*4*delta2 + 4./3.*np.exp(-spread*delta2)*delta2*spread*(-3 + 2*delta2*spread)*x**2
    else:
        res=(np.exp(-spread*(x - delta)**2)-np.exp(-spread*(x + delta)**2))*delta/(x*spread)
    return res
    
def rad_histo(data, delta, r_k, spread):
    ly=delta*0.0
    for x in data:
        ly+=r_k(x, delta, spread)
    return ly  


def get_np(fname, prefix, bsize, P, mamu, Tkelv, s, ns):   
   
    # initialises grids.
    T= Tkelv*3.1668105*10**(-6) 
    m = 1822.888* mamu
    prefix = prefix + '_'
    dq = np.zeros((bsize,3) , float)
    dqxgrid = np.linspace(-s, s, ns)
    dqygrid = np.linspace(-s, s, ns)
    dqzgrid = np.linspace(-s, s, ns)
    deltarad = dqxgrid*0.0  
    deltarad = np.sqrt(dqxgrid**2+ dqygrid**2 + dqzgrid**2)
   
    hradlist =[]
    npradlist = []

    # Defines the grid for momentum.
    pgrid = np.linspace(0.0001, np.pi /(deltarad[1]- deltarad[0]) , ns)
    pstep = np.abs(pgrid[0]-pgrid[1])
    rad_npd = pgrid*0.0
    
    #Read the end to end distances from file
    delta= np.loadtxt(fname)
    step = np.shape(delta)[0] 
   
    n_block =int(step/bsize)

    if (n_block ==0):
             print 'not enough data to build a block'
             exit()
    for x in xrange(n_block):
        dq = delta[x*bsize : (x+1)*bsize]
        dq_module = np.sqrt((dq.T[0])**2 + (dq.T[1])**2 + (dq.T[2])**2)											
     
        hrad = rad_histo(dq_module, deltarad, rad_kernel, (0.5 * T * P * m))  
        hradlist.append(hrad)
        
        # Computes the Fourier transform of the end to end vector.
        rad_npd = pgrid*0.0
        for i in range(len(pgrid)):
                 rad_npd[i] = (pgrid[i]*hrad*np.sin(pgrid[i]*deltarad)/deltarad).sum()
        npradlist.append(rad_npd) 
    
    #save the convoluted histograms of the end-to-end distances   
    avghrad = np.mean(np.asarray(hradlist), axis = 0)
    normhrad=np.sum(avghrad)
    errhrad = np.std(np.asarray(hradlist), axis = 0)/ np.sqrt(n_block)/normhrad
   
    np.savetxt(str(prefix + "rad-histo"), np.c_[deltarad, avghrad, errhrad])

   
    #save the radial n(p) and print the average value of p-square
    avgnprad = np.mean(np.asarray(npradlist), axis = 0)  
    norm=np.sum(avgnprad*pstep)
    errnprad = np.std(np.asarray(npradlist), axis = 0)/np.sqrt(n_block)/norm
    avgnprad= avgnprad/(norm)    
    np.savetxt(str(prefix + "rad-np"), np.c_[pgrid, avgnprad, errnprad])

    psqmedrad =  0.
    psqmed2rad = 0.
    for i in range(n_block):         
         psqmedrad += np.dot(pgrid**2,np.asarray(npradlist)[i,:])/np.sum(np.asarray(npradlist)[i,:])
         psqmed2rad +=  (np.dot(pgrid**2, np.asarray(npradlist)[i,:])/np.sum(np.asarray(npradlist)[i,:]))**2
    print 'av_p^2', psqmedrad/n_block, 'sigma', np.sqrt((psqmed2rad/n_block) - (psqmedrad/n_block)**2)/np.sqrt(n_block)
   

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-qfile",type=str, help="name of the end-to-end distances file")
    parser.add_argument("--prefix",type=str, default="", help="prefix for the output files")
    parser.add_argument("-bsize", type=int, default=50000, help="Specify the size of the blocks")
    parser.add_argument("-P", type=int, default= 1, help="Specify the number of beads")
    parser.add_argument("-m", type=float, default= 1.0078, help="Specify the mass of the atom in a.m.u. default is hydorgen")
    parser.add_argument("-T", type=float, default= 300, help="Specify the temperature of the system in kelvin")
    parser.add_argument("-dint", type=float, default=10, help="Specify the positive extrema of the interval to build the histogram ([-dint,dint])")
    parser.add_argument("-ns", type=float, default=1000, help="Specify the number of point to use for the histogram")
    args = parser.parse_args()

    get_np(args.qfile, args.prefix, args.bsize, args.P, args.m, args.T, args.dint, args.ns)


