#!/usr/bin/env python2

""" Instanton_interpolation.py
Reads a hessian file  and/or a positions file (xyz format) and creates an interpolation
that can be used in a further calculation

Syntax:    python  Instanton_interpolation.py  -xyz <geometry file> -h <hessian file> -n <new-beads(half-polymer)>
Example:   python  Instanton_interpolation.py  -xyz INSTANTON.xyz  -c  INSTANTON.hess -n 30

Relies on the infrastructure of i-pi, so the ipi package should
be installed in the Python module directory, or the i-pi
main directory must be added to the PYTHONPATH environment variable.

"""
#Y. Litman 2017

import os
import numpy as np
import sys
import argparse


# Y. Litman, 2017.

#You can insert the i-pi path with the following lines.
#Uncomment them and adjust the ipi_path variable

ipi_path='/home/litman/Yair/Instanton/I-PI-mc/i-pi-mc'

if not (os.path.exists(ipi_path)):
   print 'We can not find ipi in %s' %ipi_path
   print 'Please correct the path'
   sys.exit()
sys.path.insert(0, ipi_path)

from ipi.utils.io import read_file, print_file
from ipi.utils.nmtransform import nm_rescale
from ipi.utils.units import unit_to_internal 

#INPUT
parser = argparse.ArgumentParser( description="""Script for interpolate hessian and/or instanton geometry""")
parser.add_argument('-xyz','--xyz',required=True,type=str,help="Name of the instanton geometry file.")
parser.add_argument('-hess','--hessian',type=str,default='None',help="Name of the hessian file.")
parser.add_argument('-n','--nbeadsNew',required=True,default=0,  help='New number of beads (half polymer)', type=int )

args       = parser.parse_args()
input_geo  = args.xyz
input_hess = args.hessian
nbeadsNew  = args.nbeadsNew


#-----Some functions-----------------

def get_double_h(nbeads0,natoms,h0):
    """Takes nbeads, positions and hessian (only the 'physcal part') of the half polymer and 
       returns the equivalent for the full ringpolymer."""

    ii         = 3*natoms
    iii        = 3*natoms*nbeads0
    h          = np.zeros((ii , iii*2))
    h[:,0:iii] = h0

    for i in range(nbeads0):
       x = i*ii+iii
       y = ((nbeads0-1)-i)*ii
       h[:,x:x+ii] = h0[:,y:y+ii]

    return h

##################################   OPEN AND READ   ###########################################################3



if input_geo !=None:

   ipos=open(input_geo,"r")
   out=open("NEW_INSTANTON.xyz","w")
   pos=list()
   nbeads = 0
   while True:
      try:
         ret = read_file("xyz", ipos)
         pos.append(ret["atoms"])
         cell = ret["cell"]
         nbeads+=1
      except EOFError: # finished reading files
         break
   ipos.close()

   natoms=pos[0].natoms
   print ' '
   print 'We have a half ring polymer made of %i beads and %i atoms.' %(nbeads,natoms)
   print 'We will expand the ring polymer to get a half polymer of %i beads.' %nbeadsNew
   
   #Compose the half ring polymer.
   q   = np.vstack([i.q for i in pos])
   #Compose the full ring polymer.
   q2 = np.concatenate((q, np.flipud(q)), axis=0)

   #Make the rpc step
   rpc                 = nm_rescale(2*nbeads,2*nbeadsNew)
   new_q               = rpc.b1tob2(q2)[0:nbeadsNew]

   #Print
   for i in range(nbeadsNew):
     pos[0].q=new_q[i]/unit_to_internal("length", "angstrom", 1.0) #Go back to angstrom
     print_file("xyz",pos[0],cell,out,title='cell{atomic_unit}  Traj: positions{angstrom}')
     #print_file("xyz",pos[0],cell,out,title='cell  }')
   out.close()

   print 'The new Instanton geometry (half polymer) was generated'
   print 'Check NEW_INSTANTON.xyz'
   print ''

if input_hess !='None':

   print 'The new hessian is %i x %i.' %(3*natoms,natoms*3*nbeadsNew)
   if (os.path.exists(input_hess)):
       hess=open(input_hess,"r")
   else:
     print "We can't find %s " %input_hess
     sys.exit()

   out=open("NEW_HESSIAN.dat","w")
   h=np.zeros( (natoms*3)**2 * nbeads)
   aux=hess.readline().split()
   
   for i in range( (natoms*3)**2 * nbeads):
       h[i]=float(aux[i])
   h=h.reshape((natoms*3,natoms*3*nbeads))
   hess.close()

   print 'Creating matrix... '
   hessian = get_double_h(nbeads,natoms,h)

   size0 = natoms*3
   size1 = size0*(2*nbeads)
   size2 = size0*(2*nbeadsNew)

   new_h=np.zeros([size0,size2])
   for i in range(size0):
      for j in range(size0):
         h = np.array([])
         for n in range(2*nbeads):
             h = np.append(h,hessian[i,j+size0*n])
         diag=rpc.b1tob2(h)
         new_h[i,j:size2:size0] += diag

   new_h_half=new_h[:,0:size2/2]
   np.savetxt(out, new_h_half.reshape(1,new_h_half.size))

   print 'The new physical Hessian (half polymer) was generated'
   print 'Check NEW_HESSIAN.dat'
   print ''

sys.exit()
