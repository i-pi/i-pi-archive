#!/usr/bin/env python2 
"""
Computes the quantum momentum distribution of a particle given the end-to-edn distances.  
It computes both the three components of the momentum distribution and the radial function.
Moreover it computes <p^2> both in the various directions and the total contribute.
"""

import argparse
import numpy as np
import time

def kernel(x, mean=0, sigma=1):
    return np.exp(-(x- mean)**2*(0.5*sigma**2))

def histo_der(data, delta, k, mean, sigma):
    ly=delta*0.0
    ns= len(ly)
    dqstep=delta[1]-delta[0]
    for x in data:
        q = int(x/dqstep + ns/2.)
        index = q + np.arange(-int(6./dqstep/sigma), int(6./dqstep/sigma))  
        y = np.where((index < ns) & (index >= 0))   
        ly[index[y]] += -x  * k(delta[index[y]]-x, mean, sigma)           
    return ly

def histo(data, delta, k, mean, sigma):
    ly=delta*0.0
    ns= len(ly)
    dqstep=delta[1]-delta[0]
    for x in data:
        q = int(x/dqstep + ns/2.)
        index = q + np.arange(-int(6./dqstep/sigma), int(6./dqstep/sigma))  
        y = np.where((index < ns) & (index >= 0))   
        ly[index[y]] +=  k(delta[index[y]]-x, mean, sigma)           
    return ly

def rad_kernel(x, delta, spread):
   if (x <= 10**(-4)):
      res= np.exp(-spread*delta**2)*4*delta**2 + 4./3.*np.exp(-spread*delta**2)*delta**2*spread*(-3 + 2*delta**2*spread)*x**2
   else:
      res=(np.exp(-spread*(x - delta)**2)-np.exp(-spread*(x + delta)**2))*delta/(x*spread)
   return res
    
def rad_histo(data, delta, r_k, spread):
    ly=delta*0.0
    ns= len(ly)
    dqstep=delta[0]-delta[1]
    for x in data:
       q = int(x/dqstep + ns/2.)
       index = q + np.arange(-int(2000), int(2000))
       y = np.where((index < ns) & (index >= 0))  
       ly[index[y]] += r_k(x, delta[index[y]], spread)
       #ly += r_k(x, delta, spread)
    return ly  


def get_np(path, fname, bsize, P, m, Tkelv, nskip, s, ns, der):
    start = time.clock()
    print der
    # initialises grids.
    T= Tkelv*3.1668105*10**(-6) 
    dq = np.zeros((bsize,3) , float)
    dqxgrid = np.linspace(-s, s, ns)
    dqygrid = np.linspace(-s, s, ns)
    dqzgrid = np.linspace(-s, s, ns)
    deltarad = dqxgrid*0.0  
    deltarad = np.sqrt(dqxgrid**2+ dqygrid**2 + dqzgrid**2)
    dqstep= np.abs(dqxgrid[1]-dqxgrid[0])
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

    if der == False:
	    for x in xrange(n_block):
		dq = delta[x*bsize : (x+1)*bsize]
		dq_module = np.sqrt((dq.T[0])**2 + (dq.T[1])**2 + (dq.T[2])**2)											
		hx = histo(np.concatenate((dq.T[0], -dq.T[0])), dqxgrid, kernel, 0, np.sqrt(T * P * m))
		hx = ((hx + hx[::-1]) / 2.0)
		hy = histo(np.concatenate((dq.T[1], -dq.T[1])), dqygrid, kernel, 0, np.sqrt(T * P * m))
		hy = ((hy + hy[::-1]) / 2.0)
		hz = histo(np.concatenate((dq.T[2], -dq.T[2])), dqzgrid, kernel, 0, np.sqrt(T * P * m))
		hz = ((hz + hz[::-1]) / 2.0)
		hxlist.append(hx)
		hylist.append(hy)
		hzlist.append(hz)
		
		# Computes the Fourier transform of the end to end vector.
		npx = np.abs(np.fft.fftshift(np.fft.fft(hx)))
		npy = np.abs(np.fft.fftshift(np.fft.fft(hy)))
		npz = np.abs(np.fft.fftshift(np.fft.fft(hz)))
		nplistx.append(npx)
		nplisty.append(npy)
		nplistz.append(npz) 
    else:
	    for x in xrange(n_block):
		dq = delta[x*bsize : (x+1)*bsize]
		dq_module = np.sqrt((dq.T[0])**2 + (dq.T[1])**2 + (dq.T[2])**2)											
		hx = histo_der(np.concatenate((dq.T[0], -dq.T[0])), dqxgrid, kernel, 0, np.sqrt(T * P * m))
		hx = np.cumsum((hx - hx[::-1]) / 2.0)
		hy = histo_der(np.concatenate((dq.T[1], -dq.T[1])), dqygrid, kernel, 0, np.sqrt(T * P * m))
		hy = np.cumsum((hy - hy[::-1]) / 2.0)
		hz = histo_der(np.concatenate((dq.T[2], -dq.T[2])), dqzgrid, kernel, 0, np.sqrt(T * P * m))
		hz = np.cumsum((hz - hz[::-1]) / 2.0)
		hxlist.append(hx)
		hylist.append(hy)
		hzlist.append(hz)
		
		# Computes the Fourier transform of the end to end vector.
		npx = np.abs(np.fft.fftshift(np.fft.fft(hx)))
		npy = np.abs(np.fft.fftshift(np.fft.fft(hy)))
		npz = np.abs(np.fft.fftshift(np.fft.fft(hz)))
		nplistx.append(npx)
		nplisty.append(npy)
		nplistz.append(npz)
	    
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
         
    print '# number of blocks', n_block
    print '# av_px^2', psqmedx/n_block, 'sigmax', np.sqrt((psqmed2x/n_block) - (psqmedx/n_block)**2)/np.sqrt(n_block)
    print '# av_py^2', psqmedy/n_block, 'sigmay', np.sqrt((psqmed2y/n_block) - (psqmedy/n_block)**2)/np.sqrt(n_block)
    print '# av_pz^2', psqmedz/n_block, 'sigmaz', np.sqrt((psqmed2z/n_block) - (psqmedz/n_block)**2)/np.sqrt(n_block)
 
    print "# time taken (s)", time.clock()-start 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument("--path",type=str, default="", help="path of the folder conatining the end-to-end distances file")
    parser.add_argument("--fname",type=str,default="", help="name of the end-to-end distances file")
    parser.add_argument("-bsize", type=int, default=100000, help="Specify the size of the blocks")
    parser.add_argument("-P", type=int, default= 1, help="Specify the number of beads")
    parser.add_argument("-m", type=float, default= 1837, help="Specify the mass of the atom in atomic units-default is hydorgen")
    parser.add_argument("-T", type=float, default= 300, help="Specify the temperature of the system in kelvin")
    parser.add_argument("-nskip", type=int, default= 10, help="Removes the equilibration steps")
    parser.add_argument("-dint", type=float, default=15, help="Specify the positive extrema of the interval to build the histogram ([-dint,dint])")
    parser.add_argument("-ns", type=float, default=5000, help="Specify the number of point to use for the histogram")
    parser.add_argument("-der", action="store_true", default=False, help="Derives, integrates and then takes the Fourier transform")
    args = parser.parse_args()

    get_np(args.path, args.fname, args.bsize, args.P, args.m, args.T, args.nskip, args.dint, args.ns, args.der)


