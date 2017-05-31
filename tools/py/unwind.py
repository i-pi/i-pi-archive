""" trimsim.py

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
from ipi.engine.outputs import *
from ipi.engine.properties import getkey
from ipi.inputs.simulation import InputSimulation
from ipi.utils.io.inputs import io_xml
import scipy.linalg as sp
from scipy.interpolate import CubicSpline


def input_vvac(path2inputfile, mrows, stride):
    """Imports the vvac file and extracts the ."""
    try:
        dvvac=np.genfromtxt(path2inputfile, usecols=((2,3)))
        if( mrows == -1 ):
            mrows = len(dvvac)
        return dvvac[:mrows][::stride]
    except:
        print("error in inporting the vvac file",  sys.exc_info()[0])
        raise

def output_vvac(xy,path2outfile, refvvac, action):
    """Imports the vvac file and extracts the ."""
    try:
        xorg=refvvac[:,0]
        xred=xy[0]
        yred=xy[1]
        if(action == "org"):
            x=xorg
            y=CubicSpline(xred, yred, extrapolate=True)(xorg)
            np.savetxt(path2outfile,np.vstack((x, y)).T)
        elif(action == "red"):
            x=xred
            y=yred
            np.savetxt(path2outfile,np.vstack((x, y)).T)
    except:
        print("error in printing the vvac",  sys.exc_info()[0])
        raise


def Aqp(omega_0, Ap):
    """Given the free particle Ap matrix and the frequency of the harmonic oscillator, computes the full drift matrix."""
    dAqp = np.zeros(np.asarray(Ap.shape) + 1)
    dAqp[0,1] = -np.power(omega_0,2)
    dAqp[1,0] = 1
    dAqp[1:,1:] = Ap.T
    return dAqp

def Dqp(omega_0, Dp):
    """Given the free particle Dp matrix and the frequency of the harmonic oscillator, computes the full D matrix."""
    dDqp = np.zeros(np.asarray(Dp.shape) + 1)
    dDqp[1:,1:] = Dp.T
    return dDqp

def Cqp(omega_0, Ap, Dp):
    """Given the free particle Ap and Dp matrices and the frequency of the harmonic oscillator, computes the full covariance matrix."""
    dAqp = Aqp(omega_0, Ap)
    dDqp = Dqp(omega_0, Dp)
    return sp.solve_continuous_are( -dAqp, np.zeros(dAqp.shape), dDqp, np.eye(dAqp.shape[-1]))

def Cvv(omega, omega_0, Ap, Dp, dw):
    """Given the Cp and Dp matrices for a harmonic oscillator of frequency omega_0, computes the value of the Fourier transform of the velocity velocity auto-correlation function."""
    omega_0 = np.maximum(omega_0, dw)
    omega = np.maximum(omega, dw)
    dAqp = Aqp(omega_0, Ap)
    #dAqp[1,1] = np.maximum(dAqp[1,1], dAqp[1,1] +2.0*dw)
    #dAqp[1,1] += 2*dw
    dDqp = Dqp(omega_0, Dp)
    dCqp = Cqp(omega_0, Ap, Dp)
    domega2 = np.eye(dAqp.shape[-1])*np.power(omega,2)
    dAqp2 = np.dot(dAqp,dAqp)
    dCvv = np.linalg.inv(dAqp2 + domega2)
    dCvv = np.dot(np.dot(dAqp, dCvv), dCqp)
    return dCvv[1,1]/(np.pi/2.0)


def gleKernel(omega, Ap, Dp):
    """Given the Cp and Dp matrices for a harmonic oscillator of frequency omega_0, constructs the gle kernel for transformation of the velocity velocity autocorrelation function."""
    delta_omega = abs(omega[1]-omega[0])
    ngrid = len(omega) 
    dKer = np.zeros((ngrid,ngrid), float)
    for x in xrange(ngrid):
        for y in xrange(ngrid):
            dKer[x,y] = Cvv(omega[x], omega[y], Ap, Dp, delta_omega)
    return dKer*delta_omega


def ISRA(omega, ker, y, steps=500):
    """Given the thermostatted vvac spectrum and the range of frequencies, constructs the vibration density of states"""
    delta_omega = abs(omega[1]-omega[0])
    ngrid = len(omega)
    f = y
    CT = ker.T
    CTC = np.dot(ker.T, ker)
    cnvg= np.zeros((steps,2))

    for i in xrange(steps):
        f = f * np.dot(CT, y) / np.dot(CTC, f)
        ii = np.argwhere(np.isnan(f))
        f[ii] = f[ii+1]
        cnvg[i] = np.asarray((np.linalg.norm((np.dot(f,ker) - y))**2, np.linalg.norm(np.gradient(np.gradient(f)))**2))
    return f, cnvg

def main():
   
    # adds description of the program.
    parser=argparse.ArgumentParser(description="Given the parameters of a Generalized Langevin Equation and the vibrational density of states predicts the velocity-velocity autcorrelation obtained by the dynamics. Conversely, given the velocity-velocity autocorrelation function removes the disturbance affected by the thermostat and returns the underlying vibrational density of states. ")

    # adds arguments.
    parser.add_argument("-a","--action", nargs=1, choices=["conv","deconv"], default=None, help="choose conv if you want to obtain the response of the thermostat on the vibrational density of states; choose deconv if you want to obtain the micro-canonical density of states by removing the disturbance induced by the thermostat")
    parser.add_argument("-iipi", "--input_ipi", nargs=1, type=str, default=None, help="the relative path to the i-PI inputfile")
    parser.add_argument("-ivvac", "--input_vvac", nargs=1, type=str, default=None, help="the relative path to the input velocity-velocity autocorrelation function")
    parser.add_argument("-k", "--input_kernel", nargs=1, type=str, default=[None], help="the relative path to the kernel function")
    parser.add_argument("-mrows", "--maximum_rows", nargs=1, type=int, default=[-1], help="the index of the last rows to be imported from INPUT_VVAC")
    parser.add_argument("-s", "--stride", nargs=1, type=int, default=[1], help="the stride for computing the kernal")
    parser.add_argument("-ovvac", "--output_vvac", nargs=1, type=str, default=["output-vvac.data"], help="the name of the output file containing the (de)convoluted spectrum")
    parser.add_argument("-oflag","--output_flag", nargs=1, choices=["org","red"], default=["red"], help="choose orig_grid if you want OUTPUT_VVAC to have the same stride as INPUT_VVAC; choose reduced_grid if you want the OUTPUT_VVAC to have a stride of STRIDE")

    # parses arguments.
    if( len(sys.argv) > 1):
        args=parser.parse_args()
    else:
        parser.print_help()
        sys.exit()

    # stores the arguments
    path2iipi=str(args.input_ipi[-1])
    path2ivvac=str(args.input_vvac[-1])
    path2ker=args.input_kernel[-1]
    path2ovvac=str(args.output_vvac[-1])
    oflag=str(args.output_flag[-1])
    action=str(args.action[-1])
    nrows=int(args.maximum_rows[-1])
    stride=int(args.stride[-1])


    # opens & parses the input file
    ifile = open(path2iipi,"r")
    xmlrestart = io_xml.xml_parse_file(ifile) # Parses the file.
    ifile.close()

    isimul = InputSimulation()
    isimul.parse(xmlrestart.fields[0][1])
    simul = isimul.fetch()

    ttype = str(type(simul.syslist[0].motion.thermostat).__name__)

    if(ttype == "ThermoGLE"):
      Ap = simul.syslist[0].motion.thermostat.A  * 41.341373
      Cp = simul.syslist[0].motion.thermostat.C
      Dp = np.dot(Ap,Cp) + np.dot(Cp,Ap.T)

    elif(ttype == "ThermoLangevin"):
      Ap = np.asarray([1.0/simul.syslist[0].motion.thermostat.tau]).reshape((1,1)) * 41.341373
      Cp = np.asarray([1.0]).reshape((1,1))
      Dp = np.dot(Ap,Cp) + np.dot(Cp,Ap)

    ivvac=input_vvac(path2ivvac, nrows, stride)
    ix=ivvac[:,0]
    iy=ivvac[:,1]

    # computes the vvac kernel
    if (path2ker == None):
        print "# computing the kernel."
        ker = gleKernel(ix, Ap, Dp)
        np.savetxt(path2ovvac + "-ker" + str(nrows/stride)+ ".data", ker)
    else:
        print "# importing the kernel."
        ker=np.loadtxt(path2ker)

    # (de-)convolutes the spectrum
    if(action == "conv"):
        print "# printing the output spectrum."
        output_vvac((ix, np.dot(iy,ker.T)), path2ovvac, input_vvac(path2ivvac, nrows, 1), oflag)
    elif(action == "deconv"):
        print "# deconvoluting the input spectrum."
        oy, ocnvg = ISRA(ix, ker, iy)
        output_vvac((ix, oy), path2ovvac, input_vvac(path2ivvac, nrows, 1), oflag)
        np.savetxt("cnvg", ocnvg)

if __name__ == '__main__': 
   main()
