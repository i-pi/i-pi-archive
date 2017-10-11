#!/usr/bin/env python2 
description = """
Computes the quantum momentum distribution of a particle given the end-to-end distances.  
It computes both the three components of the momentum distribution and the spherically-averaged 
distribution of the proton momentum. It also computes <p^2> in the various directions, and the 
total contribute.
"""

import argparse
import numpy as np

def kernel(x, mean=0, sigma=1):
    return np.exp(-(x-mean)**2*(0.5*sigma**2))

def histo(data, delta, k, mean, sigma):
    ly=delta*0.0
    for x in data:
        ly+=k(delta-x, mean, sigma)
    return ly

def kernel3d(x, ispread2):
    return np.exp(-np.sum(x*x)*0.5*ispread2)

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


def get_np(path, fname, bsize, P, m, Tkelv, nskip, s, ns):   
   
    # initialises grids.
    T= Tkelv*3.1668105*10**(-6) 
    dq = np.zeros((bsize,3) , float)
    dqxgrid = np.linspace(-s, s, ns)
    dqygrid = np.linspace(-s, s, ns)
    dqzgrid = np.linspace(-s, s, ns)
    deltarad = dqxgrid*0.0  
    deltarad = np.sqrt(dqxgrid**2+ dqygrid**2 + dqzgrid**2)

    hxlist =[]
    hylist =[]
    hzlist =[]
    hradlist =[]

    nplistx = []
    nplisty = []
    nplistz = []
    npradlist = []

    # Defines the grid for momentum.
    pxi = -np.pi/(dqxgrid[1]-dqxgrid[0])
    pxf = +np.pi/(dqxgrid[1]-dqxgrid[0])
    pxstep = 2* np.pi/np.abs(dqxgrid[-1]-dqxgrid[0])
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

    pgrid = np.linspace(0.0001, 80, ns)
    pstep = np.abs(pgrid[0]-pgrid[1])
    rad_npd = pgrid*0.0
    
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
        dq_module = np.sqrt((dq.T[0])**2 + (dq.T[1])**2 + (dq.T[2])**2)											
     
        hx = histo(np.concatenate((dq.T[0], -dq.T[0])), dqxgrid, kernel, 0, np.sqrt(T * P * m))
        hy = histo(np.concatenate((dq.T[1], -dq.T[1])), dqygrid, kernel, 0, np.sqrt(T * P * m))
        hz = histo(np.concatenate((dq.T[2], -dq.T[2])), dqzgrid, kernel, 0, np.sqrt(T * P * m))         
        hrad = rad_histo(dq_module, deltarad, rad_kernel, (0.5 * T * P * m))  
        hxlist.append(hx)
        hylist.append(hy)
        hzlist.append(hz)
        hradlist.append(hrad)
        
        # Computes the Fourier transform of the end to end vector.
        npx = np.abs(np.fft.fftshift(np.fft.fft(hx)))
        npy = np.abs(np.fft.fftshift(np.fft.fft(hy)))
        npz = np.abs(np.fft.fftshift(np.fft.fft(hz)))
        rad_npd = pgrid*0.0
        for i in range(len(pgrid)):
             for t in range(len(deltarad)): 
                 rad_npd[i] += pgrid[i]*hrad[t]*np.sin(pgrid[i]*deltarad[t])/deltarad[t]
        nplistx.append(npx)
        nplisty.append(npy)
        nplistz.append(npz)
        npradlist.append(rad_npd) 
    
    #save the convoluted histograms of the end-to-end distances
    avghx = np.mean(np.asarray(hxlist), axis = 0)
    normhx= np.sum(avghx)
    errhx = np.std(np.asarray(hxlist), axis = 0)/ np.sqrt(n_block)/normhx
    avghy = np.mean(np.asarray(hylist), axis = 0)
    normhy= np.sum(avghy)
    errhy = np.std(np.asarray(hylist), axis = 0)/ np.sqrt(n_block)/normhy
    avghz = np.mean(np.asarray(hzlist), axis = 0)
    normhz= np.sum(avghz)
    errhz = np.std(np.asarray(hzlist), axis = 0)/ np.sqrt(n_block)/normhz
    np.savetxt(str(path + "histo.data"), np.c_[dqxgrid, avghx, errhx, dqygrid, avghy, errhy, dqzgrid, avghz, errhz])    

    avghrad = np.mean(np.asarray(hradlist), axis = 0)
    normhrad=np.sum(avghrad)
    errhrad = np.std(np.asarray(hradlist), axis = 0)/ np.sqrt(n_block)/normhrad
   
    np.savetxt(str(path + "rad-histo.data"), np.c_[deltarad, avghrad, errhrad])

    #save the resulting momentum distribution for each axes
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
    
    np.savetxt(str(path + "np.data"), np.c_[pxgrid,avgnpx/pxstep,errnpx/pxstep,avgnpy/pystep,errnpy/pystep,avgnpz/pzstep,errnpz/pzstep])   

    #print the average value of p-square for each direction
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
 

    
    #save the radial n(p) and print the average value of p-square
    avgnprad = np.mean(np.asarray(npradlist), axis = 0)  
    norm=np.sum(avgnprad*pstep)
    errnprad = np.std(np.asarray(npradlist), axis = 0)/np.sqrt(n_block)/norm
    avgnprad= avgnprad/(norm)    
    np.savetxt(str(path + "rad_np.data"), np.c_[pgrid, avgnprad, errnprad])

    psqmedrad =  0.
    psqmed2rad = 0.
    for i in range(n_block):         
         psqmedrad +=  pstep*np.dot(pgrid**2,np.asarray(npradlist)[i,:])/norm
         psqmed2rad +=  (pstep*np.dot(pgrid**2, np.asarray(npradlist)[i,:])/norm)**2
    print 'av_p^2', psqmedrad/n_block, 'sigma', np.sqrt((psqmed2rad/n_block) - (psqmedrad/n_block)**2)/np.sqrt(n_block)
   

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--path",type=str, default="", help="path of the folder conatining the end-to-end distances file")
    parser.add_argument("--fname",type=str,default="", help="name of the end-to-end distances file")
    parser.add_argument("-bsize", type=int, default=80000, help="Specify the size of the blocks")
    parser.add_argument("-P", type=int, default= 1, help="Specify the number of beads")
    parser.add_argument("-m", type=float, default= 1837, help="Specify the mass of the atom in atomic units-default is hydorgen")
    parser.add_argument("-T", type=float, default= 300, help="Specify the temperature of the system in kelvin")
    parser.add_argument("-nskip", type=int, default= 10, help="Removes the equilibration steps")
    parser.add_argument("-dint", type=float, default=10, help="Specify the positive extrema of the interval to build the histogram ([-dint,dint])")
    parser.add_argument("-ns", type=float, default=1000, help="Specify the number of point to use for the histogram")
    args = parser.parse_args()

    get_np(args.path, args.fname, args.bsize, args.P, args.m, args.T, args.nskip, args.dint, args.ns)


