#!/usr/bin/env python2

""" new_RESTART.py

Reads all the information needed from the RESTART file and
crates a new RESTART with an extended number of beads and the 
correponding extended hessian.

Syntax:
            python new_RESTART.py <name_input>  <temperature (K)> <new_beads> 
    
  Example:  python new_RESTART.py RESTART            150               64 

"""
#04Oct17

import numpy as np
import sys
import os
import math
import time

#I-PI path
ipi_path='/home/litman/Yair/Instanton/I-PI-mc/i-pi-mc'
if not (os.path.exists(ipi_path)):
   print 'We can not find ipi in %s' %ipi_path
   print 'Please correct the path'
   sys.exit()

sys.path.insert(0, ipi_path)
from ipi.engine.simulation import *
from ipi.utils.nmtransform import nm_rescale
from ipi.engine.beads import Beads

def get_doble(q0,nbeads0,natoms,u0,f0,h0):
    """Takes nbeads, positions and hessian (only the 'physcal part') of the half polymer and 
       returns the equivalent for the full ringpolymer."""
    q          = np.concatenate((q0, np.flipud(q0)), axis=0)
    u          = np.concatenate((u0, np.flipud(u0)), axis=0)
    f          = np.concatenate((f0, np.flipud(f0)), axis=0)
    nbeads     = 2*nbeads0
    ii         = 3*natoms
    iii        = 3*natoms*nbeads0

    h          = np.zeros((ii , iii*2))
    h[:,0:iii] = h0

    #diagonal block
    for i in range(nbeads0):
        x = i*ii+iii
        y = ((nbeads0-1)-i)*ii
        h[:,x:x+ii] = h0[:,y:y+ii]

    return q,h,u,f

##################################   START   ###########################################################3

#np.set_printoptions(precision=3, suppress=True, threshold=np.nan,linewidth=600)
np.set_printoptions(  threshold=np.nan)

#I/O
inputt  = sys.argv[1]
temp    = float(sys.argv[2])/315774.66
newb    = int(sys.argv[3])


print 'Start. Reading %s ...' %inputt
time0=time.time()
#
simulation = Simulation.load_from_xml(inputt, custom_verbosity='low',request_banner=False)
beads      = simulation.syslist[0].motion.beads.copy()
m          = simulation.syslist[0].motion.beads.m.copy()
nbeads     = simulation.syslist[0].motion.beads.nbeads
natoms     = simulation.syslist[0].motion.beads.natoms
hessian    = simulation.syslist[0].motion.hessian.copy()
mode       = simulation.syslist[0].motion.mode
old_u      = simulation.syslist[0].motion.old_u
old_f      = simulation.syslist[0].motion.old_f



#print hessian.shape
#print type(hessian)
#for i in range(nbeads*natoms*3):
# print i,hessian[i]

if nbeads ==1: 
   raise ValueError("We can not work from an initial geometry of only 1 bead")	
elif nbeads >newb:
   raise ValueError("This script is not prepared for contract the ring polymer")	


#Units. We use internally atomic unit.
b2a  = 0.52917721
h2K  = 3.1668152e-06
kb   = 1.0
hbar = 1.0
amu  = 1822.8885
au2K = 315774.66

#Print information
time1=time.time()
print 'End of reading. Time: %f s.'%(time1-time0) 
print '' 
print 'We have %i beads and %i atoms.' %(nbeads,natoms)
print 'We will expand the RP  to get %i beads.' %newb
print 'The new hessian is %i x %i.' %(3*natoms,newb*3*natoms)
print 'The new temperature is  %f K (%f a.u.). '%(temp*au2K,temp)
print 'Note that we are working in the %s mode.' %(mode)
print ''
print 'Creating matrix... '
    
#We create the full ring polymer and work with it. Then we take the half polymer stuff.
pos,hessian2,u,f    = get_doble(beads.q,nbeads,natoms,old_u,old_f,hessian)
hessian             = hessian2
rpc                 = nm_rescale(2*nbeads,2*newb)
new_q               = rpc.b1tob2(pos)[0:newb]
new_u               = rpc.b1tob2(u)[0:newb]
new_f               = rpc.b1tob2(f)[0:newb]


#Hessian
size0 = natoms*3
size1 = size0*(2*nbeads)
size2 = size0*(2*newb)

new_h=np.zeros([size0,size2])
for i in range(size0):
    for j in range(size0):
        h = np.array([])
        for n in range(2*nbeads):
            h = np.append(h,hessian[i,j+size0*n])
        diag=rpc.b1tob2(h)
        new_h[i,j:size2:size0] += diag

new_h_half=new_h[:,0:size2/2]
time2=time.time()
print 'We have created the new matrix in %f s.'%(time2-time1) 


##Modify simulation object
print ''
print 'Storing...'
#Hessian
simulation.syslist[0].motion.hessian=new_h_half

#Beads
newbeads       = Beads(beads.natoms, newb)
newbeads.m     = beads.m
newbeads.names = beads.names
newbeads.q     = new_q

simulation.syslist[0].beads=newbeads

#Others
simulation.syslist[0].motion.old_x=new_q
simulation.syslist[0].motion.old_u=new_u
simulation.syslist[0].motion.old_f=new_f

simulation.syslist[0].ensemble.temp=temp
simulation.step = 0
simulation.chk.store()
time3=time.time()
print 'New simulation object created in %f s.'%(time3-time2) 

print ''
print 'We are writting the new file...'
simulation.chk.write(store=False)
time4=time.time()
print 'That took  %f s.'%(time4-time3) 

print '' 
print '' 
print 'New RESTART file was created successfully.'
print 'Total time: %f s' %(time4-time0) 
print '' 


sys.exit()

