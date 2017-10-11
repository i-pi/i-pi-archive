#!/usr/bin/env python2 
description = """
Computes the momentum distribution having as input the end-to-end vectors of the open path
in atomic units. The result is the 3D distribution in atomic units, with the format
px py pz n(p)
....
"""

import argparse
import numpy as np
import time


def histo3d(data, delta, ns, cut, dqstep, mean, sigma):
    ly=delta[0,:]*0.0
    lspread = np.arange(-int(cut/sigma/dqstep), int(cut/sigma/dqstep)+1)        
    index = np.zeros(len(lspread))
    for x in data:
        qx= int(x[0]/dqstep + ns/2.)*ns*ns
        qy= int(x[1]/dqstep + ns/2.)*ns
        qz= int(x[2]/dqstep + ns/2.) 
        gx,gy,gz =np.meshgrid(qx + lspread*ns*ns, qy + lspread*ns, qz + lspread )
        index = (np.array(gx)).flatten() + (np.array(gy)).flatten() + (np.array(gz)).flatten()
        y = np.where((index < ns**3) & (index >= 0))
        ly[index[y]] += np.exp(-((delta[0, index[y]]-x[0])**2 + (delta[1, index[y]]- x[1])**2 + (delta[2, index[y]] - x[2])**2)*(0.5*sigma**2))
    return ly
    
def outer3(*vs):
    return reduce(np.multiply.outer, vs)
     
def get_np(path, fname, bsize, P, m, Tkelv, nskip, s, ns, cut):   
    start = time.clock()
    T= 3.1668105e-06*Tkelv
    dq = np.zeros((bsize,3) , float)
    dqxgrid = np.linspace(-s, s, ns)
    dqygrid = np.linspace(-s, s, ns)
    dqzgrid = np.linspace(-s, s, ns)        
    dqstep= np.abs(dqxgrid[1]-dqxgrid[0])
    halfsigma2 = 0.5*(T * P * m)
    print dqstep, 1./np.sqrt(T * P * m)
    
    data_path = str(path + fname)
    data = np.loadtxt(data_path)
    
    myh = np.zeros((ns,ns,ns))
    fx = np.zeros(ns); fy = np.zeros(ns); fz = np.zeros(ns); 
    dq = 5
    for x,y,z in data:
        qx= int(x/dqstep + ns/2.)
        qy= int(y/dqstep + ns/2.)
        qz= int(z/dqstep + ns/2.) 
        
        fx[qx-dq:qx+dq] = np.exp(-(x-dqxgrid[qx-dq:qx+dq])**2*halfsigma2)
        fy[qy-dq:qy+dq] = np.exp(-(y-dqxgrid[qy-dq:qy+dq])**2*halfsigma2)
        fz[qz-dq:qz+dq] = np.exp(-(z-dqxgrid[qz-dq:qz+dq])**2*halfsigma2)
        
        #myh[qx-dq:qx+dq,qy-dq:qy+dq,qz-dq:qz+dq] += np.einsum('i,j,k->ijk', fx[qx-dq:qx+dq], fy[qy-dq:qy+dq], fz[qz-dq:qz+dq])
        myh[qx-dq:qx+dq,qy-dq:qy+dq,qz-dq:qz+dq] += outer3(fx[qx-dq:qx+dq], fy[qy-dq:qy+dq], fz[qz-dq:qz+dq])

                
    np.savetxt(str(path + "np3d.data"), myh.flatten())    
    return
    xxx,yyy,zzz= np.meshgrid(dqxgrid, dqygrid, dqzgrid)
    qgrid= np.zeros((3,ns*ns*ns))
    qgrid[0,:] = (np.array(yyy)).flatten()
    qgrid[1,:] = (np.array(xxx)).flatten()
    qgrid[2,:] = (np.array(zzz)).flatten()

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
 
    px,py,pz= np.meshgrid(pxgrid, pygrid, pzgrid)
    pgrid= np.zeros((3,ns*ns*ns))
    pgrid[0,:] = (np.array(py)).flatten()
    pgrid[1,:] = (np.array(px)).flatten()
    pgrid[2,:] = (np.array(pz)).flatten() 
    dpstep= np.abs(pxgrid[1]-pxgrid[0])

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
        h3d = histo3d(np.concatenate((dq, -dq)), qgrid, ns, cut, dqstep, 0, np.sqrt(T * P * m)) 
        # Computes the 3D Fourier transform of the histogram of the end-to-end distances.
        npvec3d= np.abs(np.fft.fftshift(np.fft.fftn(np.reshape(h3d, (ns, ns, ns)))))       
        nplist3d.append(npvec3d.flatten())
        
    avgnp3d= np.mean(np.asarray(nplist3d), axis = 0)
    norm= sum(avgnp3d*dpstep**3)
    avgnp3d= avgnp3d/norm
    np.savetxt(str(path + "np3d.data"), np.c_[pgrid.T, avgnp3d])
    print time.clock()-start 

if __name__ == '__main__':
    # get_npvec delta-file-name --prefix [output-prefix]
    # the outputs will be like output-prefix.np3d, output-prefix.delta3d 
    # -dint should be selected automatically if not given as maximum value seen in the input *1.1
    # -ns default should be like 0, if detects it's zero, automatically set so that the grid spacing is sigma of the kernel
    # the threshold should be by default like 6 sigma
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--path",type=str, default="", help="path of the folder conatining the end-to-end distances file")
    parser.add_argument("--fname",type=str,default="", help="name of the end-to-end distances file")
    parser.add_argument("-bsize", type=int, default=50000, help="Specify the size of the blocks")
    parser.add_argument("-P", type=int, default= 1 , help="Specify the number of beads")
    parser.add_argument("-m", type=float, default= 1837 , help="Specify the mass of the atom in atomic units-default is hydorgen")
    parser.add_argument("-T", type=float, default= 300 , help="Specify the temperature of the system in kelvin")
    parser.add_argument("-nskip", type=int, default= 500 , help="Removes the equilibration steps")
    parser.add_argument("-dint", type=float, default=10.0, help="Specify the positive extrema of the interval to build the histogram ([-dint,dint])")
    parser.add_argument("-ns", type=int, default=150, help="Specify the number of point to use for the histogram")
    parser.add_argument("-cut", type=int, default=15, help="Specify the size of the grid around a specific point in units of sigma")
    args = parser.parse_args()

    get_np(args.path, args.fname, args.bsize, args.P, args.m, args.T, args.nskip, args.dint, args.ns, args.cut)


