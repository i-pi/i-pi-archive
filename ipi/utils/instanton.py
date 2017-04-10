"""Contains simple helper algorithms for minimization.

Copyright (C) 2013, Joshua More and Michele Ceriotti

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http.//www.gnu.org/licenses/>.

Algorithms implemented by Yair Litman and Mariana Rossi, 2017

Functions:
        TODO

"""

__all__ = ["Instanton", ]

import numpy as np
from ipi.utils.messages import verbosity, info
from ipi.utils import units
from ipi.utils.depend import *
from ipi.utils.mintools import nichols, Powell



def Instanton(x0, f0,f1, h0, update,asr, im,gm, big_step):
    """ Input: x0 = previous positions
               f0 = previous physical forces
               f1 = previous spring forces
               h0 = previous hessian
           update = how to update the hessian
               im = instanton mapper
               gm = gradiente mapper
         big_step = limit on step length"""


    # Project out rotations and translation from the Hessian
    # Find new movement direction

    # Note that the dynmax.size < h0
    d, dynmax = clean_hessian(h0, im.dbeads.q, im.dbeads.natoms, im.dbeads.nbeads, im.dbeads.m, im.dbeads.m3,asr)
    d_x = nichols(f0,f1, d,dynmax,im.dbeads.m3, big_step)

    # Rescale step
    d_x_norm = np.linalg.norm(d_x)
    info(" @Instanton: Step norm = %g" % d_x_norm, verbosity.low)
    if d_x_norm > big_step:
        info(" @Instanton: Attempted step norm = %g, scaled down to %g" % (d_x_norm, big_step), verbosity.low)
        d_x *= big_step / d_x_norm


    # Make movement  and get new energy (u)  and forces(f) using mapper
    x = x0 + d_x
    u, g1 = im(x)
    u, g2 = gm(x)
    f = -(g1+g2)

    # Update hessian
    if update == 'powell':
        d_g = np.subtract((f0+f1), f) #Gradient dif
        Powell(d_x.flatten(), d_g.flatten(), h0)
    elif update == 'true':
        get_hessian(h0,gm,x)
        if gm.dbeads.nbeads != 1:
            add_spring_hessian(im,h0)
        u, g = gm(x) #To update pos and for values in the mapper. (Yes, it is stricly unnecesary but keeps the code simple)

#-------------------------------------------------------------------------------------------------------------
def get_hessian(h,gm,x0,d=0.01):

    info(" @GEOP: Computing hessian" ,verbosity.low)
    ii = gm.dbeads.natoms * 3
    h[:]=np.zeros((h.shape),float)


    for j in range(ii):
        x = x0.copy()

        x[:, j] = x0[:, j] + d
        e, f1 = gm(x)

        x[:, j] = x0[:, j] - d
        e, f2 = gm(x)
        g = (f1 - f2) / (2 * d)

        for i in range(gm.dbeads.nbeads):
            h[j + i * ii, i * ii:(i + 1) * ii] = g[i, :]


def add_spring_hessian(im,h):
    """ Add spring terms to the extended hessian
        """
    if im.dbeads.nbeads == 1:
        return

    ii = im.dbeads.natoms * 3

    # Spring part
    h_sp = im.dbeads.m3[0] * im.omega2
    # Diagonal
    diag = np.diag(2.0 * h_sp)

    for i in range(0, im.dbeads.nbeads):
        h[i * ii:(i + 1) * ii, i * ii:(i + 1) * ii] += diag

    # Non-Diagonal
    ndiag = np.diag(-h_sp)
    #quasi-band
    for i in range(0, im.dbeads.nbeads - 1):
        h[i * ii:(i + 1) * ii, (i + 1) * ii:(i + 2) * ii] += ndiag
        h[(i + 1) * ii:(i + 2) * ii, i * ii:(i + 1) * ii] += ndiag

    # Corner
    if im.mode=='full':
        h[ 0 :ii , (im.dbeads.nbeads-1)*ii:(im.dbeads.nbeads) * ii   ] += ndiag
        h[(im.dbeads.nbeads - 1) * ii:(im.dbeads.nbeads) * ii, 0:ii  ] += ndiag

def expand_hessian(h0,h,natoms):
    """ Expand initial hessian. Copy the orginal and put it in the diagonal blocks.
           IN    h0      = initial hessian
                 h       = # beads
                 natoms = # atoms
    """
    info(" @GEOP: We are in expand_hessian", verbosity.low)

    stop
    #ii = 3*natoms
    #hessian = np.zeros( (n*natoms*3,n*natoms*3) )

    #for i in range (n):
    #    hessian[ i*ii:(i+1)*ii , i*ii:(i+1)*ii  ] = h

    #return hessian

def get_imvector(h,  m3):
    """ Compute eigenvector  corresponding to the imaginary mode
            IN     h      = hessian
                   m3     = mass vector (dimension = 1 x 3*natoms)
                   nbeads = nbeads
            INTERNAL
                   mm = "mass weighted matrix"
                   hm = mass weighted hessian matrix
                   d  = mass weighted hessian matrix eigenvalues
                   w  = mass weighted hessian matrix eigenvector (in columns)
              TODO
        """
    ii = m3.size
    if h.size != m3.size**2:
        raise ValueError("@Get_imvector. Inital hessian size does not match system size.")
    m   = 1.0 / (m3 ** 0.5)
    mm  = np.outer(m, m)
    hm  = np.multiply(h, mm)

    #Simmetrize to use linalg.eigh
    hmT = hm.T
    hm  = (hmT+hm)/2.0

    d, w = np.linalg.eigh(hm)
    freq = np.sign(d) * np.absolute(d) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17)

    info(" @GEOP: 1 frequency %4.1f cm^-1" % freq[0], verbosity.low)
    info(" @GEOP: 2 frequency %4.1f cm^-1" % freq[1], verbosity.low)
    info(" @GEOP: 3 frequency %4.1f cm^-1" % freq[2], verbosity.low)
    if freq[0] > -80 and freq[0] < 0:
        info(" @GEOP: Warning, small negative frequency %4.1f cm^-1" % freq, verbosity.low)
    elif freq[0] > 0:
        raise ValueError("@GEOP: The smallest frequency is positive. We aren't in a TS. Please check your hessian")

    imv = w[:, 0] * (m3[:] ** 0.5)
    imv = imv/np.linalg.norm(imv)
    return imv

#-------------------------------------------------------------------------------------------------------------------
#Adapted from ipi/engine/motion/phonons.py apply_asr

def clean_hessian(h,q, natoms,nbeads,m,m3,asr):
    """
        Removes the translations and rotations modes.
        IN  h      = hessian
            q      =
            natoms =
            nbeads =
            m      =
            m3     =
        OUT d      = non zero eigenvalues of the dynmatrix
            w      = dynmatrix with the external modes projected out"""

    #Set some util things
    ii = natoms * nbeads
    mm = np.zeros((nbeads, natoms))
    for i in range(nbeads):
        mm[i] = m
    mm     = mm.reshape(ii)
    ism    = m3.reshape(ii * 3) ** (-0.5)
    ismm   = np.outer(ism, ism)
    dynmat = np.multiply(h, ismm)

    if asr =='none':
        hm=dynmat
    else:
       #Computes the centre of mass.
       com = np.dot(np.transpose(q.reshape((ii, 3))), mm) / mm.sum()
       qminuscom = q.reshape((ii, 3)) - com

       if asr =='poly':
           #Computes the moment of inertia tensor.
           moi  = np.zeros((3,3), float)
           for k in range(ii):
               moi-=np.dot(np.cross(qminuscom[k],np.identity(3)),np.cross(qminuscom[k],np.identity(3)))*mm[k]

           U=(np.linalg.eig(moi))[1] #eigenvector in columns
           R=np.dot(qminuscom,U)
           D=np.zeros((6,3*ii),float)

           #Computes the vectors along translations and rotations.
           #Translations
           D[0]=np.tile([1,0,0],ii)/ism
           D[1]=np.tile([0,1,0],ii)/ism
           D[2]=np.tile([0,0,1],ii)/ism
           #Rotations
           for i in range(3*ii):
               iatom=i/3
               idof=np.mod(i,3)
               D[3,i]=(R[iatom,1]*U[idof,2] - R[iatom,2]*U[idof,1])/ism[i]
               D[4,i]=(R[iatom,2]*U[idof,0] - R[iatom,0]*U[idof,2])/ism[i]
               D[5,i]=(R[iatom,0]*U[idof,1] - R[iatom,1]*U[idof,0])/ism[i]


           for k in range(6):
               D[k] = D[k] / np.linalg.norm(D[k])
            #Computes the transformation matrix.
           transfmatrix = np.eye(3*ii) - np.dot(D.T,D)
           hm = np.dot(transfmatrix.T,np.dot(dynmat,transfmatrix))

       elif asr == 'crystal':
           # Computes the vectors along translations.
           # Translations
           D = np.zeros((6, 3 * ii), float)
           D[0] = np.tile([1, 0, 0], ii) / ism
           D[1] = np.tile([0, 1, 0], ii) / ism
           D[2] = np.tile([0, 0, 1], ii) / ism

           for k in range(3):
               D[k] = D[k] / np.linalg.norm(D[k])
           # Computes the transformation matrix.
           transfmatrix = np.eye(3 * ii) - np.dot(D.T, D)
           hm = np.dot(transfmatrix.T, np.dot(dynmat, transfmatrix))

    ##Simmetrize to use linalg.eigh
    hmT = hm.T
    hm  = (hmT+hm)/2.0
    d, w = np.linalg.eigh(hm)

    #Count
    dd = np.sign(d) * np.absolute(d) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17) #convert to cm^-1

    #Zeros
    condition = np.abs(dd) < 0.01 #Note that dd[] units are cm^1
    nzero= np.extract(condition, dd)

    if nzero.size != 6:
         info(" @GEOP: Warning, we have %d 'zero' frequencies" %nzero.size, verbosity.low)

    #Negatives
    condition = dd < -0.01 #Note that dd[] units are cm^1
    nneg = np.extract(condition, dd)
    info(" @GEOP: We have %d 'neg' frequencies " % (nneg.size), verbosity.low)

    # Now eliminate external degrees of freedom from the dynmatrix
    d = np.delete(d,range(nneg.size,nneg.size+nzero.size))
    w = np.delete(w, range(nneg.size,nneg.size+nzero.size),axis=1)

    info(" @Clean hessian", verbosity.low)
    #print np.sign(d) * np.absolute(d).T ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17) #convert to cm^-1
    #for i in range(mp.dbeads.q.size):
    #    np.set_printoptions(precision=4, suppress=True)
    #    print h[i, i]
    return d,w

def print_instanton(h,gm,im,asr,rates):
    # Compute action
    action1 = gm.dforces.pot / (im.temp * im.dbeads.nbeads * units.Constants.kb)
    action2 = im.pot / (im.temp * im.dbeads.nbeads * units.Constants.kb)
    # Note that for the half polymer the factor 2 cancels out (One factor in *.pot and one factor in *.nbeads)
    print  'ACTION1', action1 / units.Constants.hbar
    print  'ACTION2', action2 / units.Constants.hbar
    print  'ACTION', (action1 + action2) / units.Constants.hbar

    if im.mode =='full':
        d,w = clean_hessian(h,im.dbeads.q,im.dbeads.natoms,im.dbeads.nbeads,im.dbeads.m,im.dbeads.m3,asr)
        np.set_printoptions(precision=4, suppress=True)
        print "Final  lowest six frequencies"
        print np.sign(d[0:6]) * np.absolute(d[0:6]) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17)  # convert to cm^-1
    elif im.mode =='half':
        im.mode ='full'
        get_hessian(h, gm, im.dbeads.q)
        add_spring_hessian(im, h)
        hbig=get_doble_hessian(h,im)
        q,nbeads,m,m3 = get_doble(im.dbeads.q, im.dbeads.nbeads, im.dbeads.m, im.dbeads.m3)
        d, w = clean_hessian(hbig, q, im.dbeads.natoms, nbeads, m, m3, asr)
        np.set_printoptions(precision=4, suppress=True)
        print "Final  lowest six frequencies"
        print np.sign(d[0:6]) * np.absolute(d[0:6]) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17)  # convert to cm^-1


    if rates == 'true':
        rates(im,gm)

def get_doble_hessian(h0,im):
    nbeads = im.dbeads.nbeads
    natoms = im.dbeads.natoms
    ii     = 3*natoms
    iii    = 3*natoms*nbeads

    #np.set_printoptions(precision=4, suppress=True)

    h = np.zeros((iii*2, iii*2))
    h[0:iii,0:iii]=h0


    #diagonal block
    for i in range(nbeads):
        x = i*ii+iii
        y = ((nbeads-1)-i)*ii
        h[x:x+ii,x:x+ii] = h0[y:y+ii,y:y+ii]

    #off-diagonal block-manual copy
    #for i in range(nbeads-1):
    #    x  = i * ii + iii
    #    xx = (i+1) * ii + iii
    #    y  = ((nbeads - 1) - (i+1)) * ii
    #    yy = ((nbeads - 1) - i) * ii
    #    h[x:x + ii, xx:xx + ii] = h0[y:y + ii, yy:yy + ii]
    #    h[xx:xx + ii, x:x + ii] = h0[yy:yy + ii, y:y + ii]
    #    print x, xx, y, yy

    # Spring part.
    h_sp = im.dbeads.m3[0] * im.omega2

    ndiag = np.diag(-h_sp)
    # quasi-band

    for i in range(0, 2*im.dbeads.nbeads - 1):
        h[i * ii:(i + 1) * ii, (i + 1) * ii:(i + 2) * ii] += ndiag
        h[(i + 1) * ii:(i + 2) * ii, i * ii:(i + 1) * ii] += ndiag

    # Corner
    h[0:ii, (2*nbeads - 1) * ii:(2*nbeads) * ii] += ndiag
    h[(2*nbeads - 1) * ii:(2*nbeads) * ii, 0:ii] += ndiag

    return h

def get_doble(q0,nbeads,m,m3):

    q=np.concatenate((q0, np.flipud(q0)), axis=0)
    m3=np.concatenate((m3, m3), axis=0)
    return q,2*nbeads,m,m3

def rates(im,gm):
    info(" @GEOP: We are going to do the final post " , verbosity.low)
    pass
