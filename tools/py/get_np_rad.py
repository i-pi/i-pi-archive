#!/usr/bin/env python2 
description = """
Computes the quantum momentum distribution of a particle given the end-to-end distances.  
It computes both the three components of the momentum distribution and the spherically-averaged 
distribution of the proton momentum. It also computes <p^2> in the various directions, and the 
total contribute.
"""

import argparse
import numpy as np

def rad_kernel(x, delta, is2half):
    # computes the kernel that describes the contribution to the radial average n(Delta) for a given spread factor 1/(2sigma^2)
    if (x <= 1.0e-4): # use a Taylor expansion to get a stable result for small x
        delta2 = delta**2
        res= np.exp(-is2half*delta2)*4*delta2 + 4./3.*np.exp(-is2half*delta2)*delta2*is2half*(-3 + 2*delta2*is2half)*x**2
    else:
        res=(np.exp(-is2half*(x - delta)**2)-np.exp(-is2half*(x + delta)**2))*delta/(x*is2half)
    return res

def rad_f_kernel(x, delta, is2half):
    # computes the kernel that describes the contribution to the radial average n(Delta) for a given spread factor 1/ispread
    if (x <= 1.0e-4): # use a Taylor expansion to get a stable result for small x
        res=  delta**2*(4*delta*is2half*x*(5 + is2half*(-5 + 2*delta**2*is2half)*x**2))/(15.*np.exp(delta**2*is2half))
    else:
        res = ((1 + 2*delta*is2half*x)*np.exp(-is2half*(delta + x)**2) + np.exp(4*delta*is2half*x-is2half*(delta + x)**2)*(-1 + 2*delta*is2half*x))/(4.*is2half**2*x**2)
    return res
    
def rad_histo(data, delta, r_k, spread):
    ly=delta*0.0
    for x in data:
        ly+=r_k(x, delta, spread)
    return ly  

def rad_fhisto(data, fdata, delta, r_k, spread):
    ly=delta*0.0
    for i in xrange(len(data)):
        ly+=r_k(data[i], delta, spread)*fdata[i]/data[i]
    return ly  
    
def get_np(fname, ffile, prefix, bsize, P, mamu, Tkelv, s, ns, skip):   
   
    # initialises grids.
    T= Tkelv*3.1668105*10**(-6) 
    m = 1822.888* mamu    
    prefix = prefix + '_'
    
    #Read the end to end distances from file
    delta= np.loadtxt(fname)[skip:]
    step = np.shape(delta)[0] 
    if ffile!="": 
        fdelta = np.loadtxt(ffile)[skip:]
    else:
        fdelta = None    
    
    if s<= 0 : 
        s = np.sqrt(np.max(np.sum(delta**2,axis=1)))*4    
    if ns<=0:
        ns = int(s*np.sqrt(T * P * m))*2+1
    deltarad = np.linspace(0, s, ns)    
   
    hradlist = []
    dhradlist = []
    npradlist = []
    p2list = []

    # Defines the grid for momentum.
    pgrid = np.linspace(0, np.pi /(deltarad[1]- deltarad[0]) , ns)
    pstep = np.abs(pgrid[0]-pgrid[1])
    rad_npd = pgrid*0.0    
    
    n_block =int(step/bsize)

    if (n_block ==0):
             print 'not enough data to build a block'
             exit()
    for x in xrange(n_block):
        dq = delta[x*bsize : (x+1)*bsize]
        dq_module = np.sqrt(np.sum(dq*dq,axis=1))
     
        hrad = rad_histo(dq_module, deltarad, rad_kernel, (0.5 * T * P * m))  
        hradlist.append(hrad)
        
        if not fdelta is None:
            df = fdelta[x*bsize : (x+1)*bsize]
            dfq = np.sum(df*dq,axis=1)
            dhrad = rad_fhisto(dq_module, dfq, deltarad, rad_f_kernel, (0.5 * T * P * m))  
            dhradlist.append(dhrad)
        
        # Computes the Fourier transform of the end to end vector.
        rad_npd = pgrid*0.0
        for i in range(len(pgrid)):
            # computes the FT by hand, and also takes care of the zero limit
            fsin = np.sin(pgrid[i]*deltarad);
            fsin[1:] /= deltarad[1:]
            fsin[0] = pgrid[i]
            rad_npd[i] = (pgrid[i]*hrad*fsin).sum()
        npradlist.append(rad_npd)

        #Computes the second moment. 
        p2list.append(np.dot(pgrid**2, rad_npd))
    
    #saves the convoluted histograms of the end-to-end distances   
    avghrad = np.sum(np.asarray(hradlist), axis = 0)
    normhrad = np.sum(avghrad)
    errhrad = np.std(np.asarray(hradlist), axis = 0) * np.sqrt(n_block)
   
    np.savetxt(str(prefix + "rad-histo"), np.c_[deltarad, avghrad / normhrad, errhrad / normhrad])

    if not fdelta is None:
        avgdhrad = np.mean(np.asarray(dhradlist), axis = 0)
        errdhrad = np.std(np.asarray(dhradlist), axis = 0)/ np.sqrt(n_block)
        np.savetxt(str(prefix + "rad-dhisto"), np.c_[deltarad, avgdhrad, errdhrad])
   
    #saves the radial n(p) and print the average value of p-square
    avgnprad = np.sum(np.asarray(npradlist), axis = 0)  
    normnprad = np.sum(avgnprad)
    errnprad = np.std(np.asarray(npradlist), axis = 0) * np.sqrt(n_block)
    np.savetxt(str(prefix + "rad-np"), np.c_[pgrid, avgnprad / normnprad, errnprad / normnprad])

    print 'av_p^2', np.sum(np.asarray(p2list) / normnprad), '+/-', np.std(np.asarray(p2list) / normnprad) * np.sqrt(n_block)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-qfile",type=str, help="name of the end-to-end distances file")
    parser.add_argument("-ffile",type=str, default="", help="name of the forces file")
    parser.add_argument("--prefix",type=str, default="out", help="prefix for the output files")
    parser.add_argument("-bsize", type=int, default=50000, help="Specify the size of the blocks")
    parser.add_argument("-P", type=int, default= 1, help="Specify the number of beads")
    parser.add_argument("-m", type=float, default= 1.0078, help="Specify the mass of the atom in a.m.u. default is hydorgen")
    parser.add_argument("-T", type=float, default= 300, help="Specify the temperature of the system in kelvin")
    parser.add_argument("-dint", type=float, default=0, help="Specify the positive extrema of the interval to build the histogram ([-dint,dint])")
    parser.add_argument("-ns", type=float, default=0, help="Specify the number of point to use for the histogram")
    parser.add_argument("-skip", type=int, default=0, help="Specify the number of points to be skipped")
    
    args = parser.parse_args()

    get_np(args.qfile, args.ffile, args.prefix, args.bsize, args.P, args.m, args.T, args.dint, args.ns, args.skip)


