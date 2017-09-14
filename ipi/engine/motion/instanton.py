"""
Contains classes for different geometry optimization algorithms.

TODO

Algorithms implemented by Michele Ceriotti and Benjamin Helfrecht, 2015
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np
from scipy import linalg
import time

from ipi.engine.motion import Motion
from ipi.utils.depend import depstrip, dobject
#from ipi.utils.depend import *
from ipi.utils.softexit import softexit
from ipi.utils.messages import verbosity, info
from ipi.utils import units
from ipi.utils.mintools import nichols, Powell

__all__ = ['InstantonMotion']

class InstantonMotion(Motion):
    """Instanton motion class.

    Attributes:
        mode: minimization algorithm to use
        biggest_step: max allowed step size 
        old_force: force on previous step
        hessian: stored  Hessian matrix 
        tolerances:
            energy: change in energy tolerance for ending minimization
            force: force/change in force tolerance foe ending minimization
            position: change in position tolerance for ending minimization}
        delta:
        hessian_init:
        hessian_asr:
        final_rates:
    """

    def __init__(self, fixcom=False, fixatoms=None,
                 mode="instanton",
                 biggest_step=0.5,
                 old_pos=np.zeros(0, float),
                 old_pot= np.zeros(0, float),
                 old_force=np.zeros(0, float),
                 hessian=np.eye(0, 0, 0, float),
                 tolerances={"energy": 1e-6, "force": 1e-6, "position": 1e-6},
                 delta=np.zeros(0, float),
                 hessian_init=None,
                 hessian_update=None,
                 hessian_asr=None,
                 final_rates=None):

        """Initialises InstantonMotion.
        """

        super(InstantonMotion, self).__init__(fixcom=fixcom, fixatoms=fixatoms)

        # Optimization Options
        self.mode           = mode
        self.big_step       = biggest_step
        self.tolerances     = tolerances

        self.old_x          = old_pos
        self.old_u          = old_pot
        self.old_f          = old_force

        #
        self.optimizer      = InstantonOptimizer()
        self.hessian        = hessian
        self.hessian_init   = hessian_init
        self.hessian_update = hessian_update
        self.hessian_asr    = hessian_asr
        self.final_rates    = final_rates
        self.delta          = delta


    def bind(self, ens, beads, nm, cell, bforce, prng):
        """Binds beads, cell, bforce and prng to InstantonMotion

            Args:
            beads: The beads object from whcih the bead positions are taken.
            nm: A normal modes object used to do the normal modes transformation.
            cell: The cell object from which the system box is taken.
            bforce: The forcefield object from which the force and virial are taken.
            prng: The random number generator object which controls random number generation.
        """

        super(InstantonMotion,self).bind(ens, beads, nm, cell, bforce, prng)
        # Binds optimizer
        self.optimizer.bind(self)

    def step(self, step=None):
        self.optimizer.step(step)

class GradientMapper(object):

    """Creation of the multi-dimensional function that will be minimized.
    Used in the BFGS and L-BFGS minimizers.

    Attributes:
        dbeds:   copy of the bead object
        dcell:   copy of the cell object
        dforces: copy of the forces object
    """

    def __init__(self):
        pass

    def bind(self, dumop):
        self.dbeads = dumop.beads.copy()
        self.dcell = dumop.cell.copy()
        self.dforces = dumop.forces.copy(self.dbeads, self.dcell)

    def set_pos(self,x):
        """Set the positions """
        self.dbeads.q = x

    def __call__(self,x):
        """computes energy and gradient for optimization step"""

        self.dbeads.q = x
        e = self.dforces.pot   # Energy
        g = -self.dforces.f   # Gradient

        return e, g


class InstantonMapper(object):
    """Creation of the multi-dimensional function to compute full or half ring polimer.

    Attributes:
        dbeds:   copy of the bead object
    """

    def __init__(self):
        self.pot = None
        self.f   = None
        pass

    def bind(self, dumop):
        self.dbeads = dumop.beads.copy()
        self.temp   = dumop.temp
        self.mode   = dumop.mode
        if self.mode=='full':
            self.omega2 = (self.temp * self.dbeads.nbeads * units.Constants.kb / units.Constants.hbar) ** 2
        elif self.mode =='half':
            self.omega2 = (self.temp * (2*self.dbeads.nbeads) * units.Constants.kb / units.Constants.hbar) ** 2

    def save(self,e,g):
        self.pot = e
        self.f   = -g

    def __call__(self, x):
        """Computes spring energy and gradient for instanton optimization step"""
        if x.shape[0]==1: #only one beads
           e=0.0
           g=np.zeros(x.shape[1])
           self.save(e, g)

        if self.mode=='full':
            self.dbeads.q = x.copy()
            e=self.dbeads.vpath*self.omega2
            g=-self.dbeads.fpath*self.omega2
            self.save(e,g)

        elif self.mode=='half':
            self.dbeads.q = x.copy()
            e = 0.00
            g = np.zeros(self.dbeads.q.shape, float)
            g2 = np.zeros(self.dbeads.q.shape, float)
            for i in range(self.dbeads.nbeads - 1):
                #e += (np.sum( self.dbeads.m3[0] * self.omega2 * 0.5 *(self.dbeads.q[i + 1, :] - self.dbeads.q[i, :] ) ** 2))
                dq =self.dbeads.q[i + 1, :] - self.dbeads.q[i, :]
                e += self.omega2 * 0.5 *np.dot(self.dbeads.m3[0]*dq,dq)
            for i in range(0, self.dbeads.nbeads - 1):
                g[i, :] += self.dbeads.m3[i,:] * self.omega2 * (self.dbeads.q[i, :] - self.dbeads.q[i + 1, :])
                g2[i, :] += self.dbeads.m3[i, :] * (self.dbeads.q[i, :] - self.dbeads.q[i + 1, :])
            for i in range(1, self.dbeads.nbeads ):
                g[i, :] += self.dbeads.m3[i,:] * self.omega2 * (self.dbeads.q[i, :] - self.dbeads.q[i - 1, :])
                g2[i, :] += self.dbeads.m3[i, :] * (self.dbeads.q[i, :] - self.dbeads.q[i - 1, :])

            self.save(e,g)

        return e, g

class InstantonOptimizer(dobject):
    """ INSTANTON """

    def __init__(self):

        """initialises object for LineMapper (1-d function) and for GradientMapper (multi-dimensional function) """

        self.gm           = GradientMapper()
        self.im           = InstantonMapper()

    def bind(self, geop):

        """
        bind optimization options and call bind function of Mappers (get beads, cell,forces)
        check wheter force size,  Hessian size from previous step match system size
        """
        if geop.beads.nbeads == 1 and geop.mode=='half':
            raise ValueError("Half mode is not compatible with nbeads = 1")

        self.mode       = geop.mode
        self.beads      = geop.beads
        self.cell       = geop.cell
        self.forces     = geop.forces
        self.fixcom     = geop.fixcom
        self.fixatoms   = geop.fixatoms


        #The resize action must be done before the bind
        if geop.old_x.size != self.beads.q.size:
            if geop.old_x.size == 0:
                geop.old_x =np.zeros((self.beads.nbeads,3*self.beads.natoms), float)
            else:
                raise ValueError("Old positions size does not match system size")
        if geop.old_u.size != 1:
            if geop.old_u.size == 0:
                geop.old_u = np.zeros(1, float)
            else:
                raise ValueError("Old potential energy size does not match system size")
        if geop.old_f.size != self.beads.q.size:
            if geop.old_f.size == 0:
                geop.old_f = np.zeros((self.beads.nbeads,3*self.beads.natoms), float)
            else:
                raise ValueError("Old forces size does not match system size")

        self.old_x      = geop.old_x
        self.old_u      = geop.old_u
        self.old_f      = geop.old_f
        self.temp       = geop.ensemble.temp


        if geop.ensemble.temp == -1.0 or geop.ensemble.temp == 1.0: #This is due to a little inconsistency on the default value
            if self.beads.nbeads != 1:
                raise ValueError("Temperature must be specified for an Instanton calculation ")

        self.tolerances     = geop.tolerances
        self.big_step       = geop.big_step
        self.delta          = geop.delta


        self.hessian_update  = geop.hessian_update
        self.hessian_asr     = geop.hessian_asr
        self.hessian_init    = geop.hessian_init
        self.final_rates     = geop.final_rates

        self.gm.bind(self)
        self.im.bind(self)

        #Hessian
        self.initial_hessian = None

        if geop.hessian.size != (self.beads.q.size * self.beads.q.size):
            if geop.hessian.size == 0: #Hessian not provided
                if self.beads.nbeads == 1:
                    if geop.hessian_init == 'true':
                        info(" Initial classical hessian is not provided. We are going to compute it.", verbosity.low)
                        geop.hessian = np.zeros((self.beads.q.size, self.beads.q.size))
                    else:
                        #If nbeads>1 you have to provide a initial guess for the hessian
                        raise ValueError("Hessian_init != 'True'.  An initial hessian must be provided")
                else:
                    raise ValueError("nbeads >1. An initial hessian must be provided")
            elif geop.hessian.size <=  (self.beads.q.size * self.beads.q.size):
                self.initial_hessian = geop.hessian.copy()
                geop.hessian         = np.zeros((self.beads.q.size, self.beads.q.size),float)
            else:
                raise ValueError("Hessian size does not match system size.")

        self.hessian = geop.hessian


    def step(self, step=None):
        """ Does one simulation time step.
            Attributes:
            qtime: The time taken in updating the positions.
        """
        self.qtime = -time.time()
        info("\nMD STEP %d" % step, verbosity.medium)#ALBERTO
        #np.set_printoptions(precision=10, suppress=True, threshold='nan')


        if step == 0:
            info(" @GEOP: Initializing INSTANTON", verbosity.low)

            if self.beads.nbeads == 1:
                info(" @GEOP: Classical TS search", verbosity.low)
                if self.hessian_init == 'true':
                    get_hessian(self.hessian, self.gm, self.beads.q)
                u, g = self.im(self.beads.q) #Init im
            else:
                if ((self.beads.q - self.beads.q[0]) == 0).all():  # If the coordinates in all the imaginary time slices are the same
                    info(" @GEOP: We stretch the initial geometry with an 'amplitud' of %4.2f" % self.delta, verbosity.low)
                    imvector = get_imvector(self.initial_hessian,self.beads.m3[0].flatten())
                    for i in range(self.beads.nbeads):
                        if self.mode == 'full':
                            self.beads.q[i, :] += self.delta * np.cos((i+1) *2.0*np.pi / float(self.beads.nbeads)) * imvector[:]
                        elif self.mode== 'half':
                           self.beads.q[i, :] += self.delta * np.cos(i * np.pi / float(self.beads.nbeads-1)) * imvector[:]
                           # print 'ALBERTO-line'
                           # n=self.beads.nbeads
                           # self.beads.q[:,0] =np.linspace(-1,1,n)


                    if self.hessian_init !='true':
                        info(" @GEOP: Hessian_init isn't true but we have stretched the polymer so we are going to compute the initial hessian anyway", verbosity.low)
                    self.hessian_init ='true'
                else:
                    info(" @GEOP: Starting from the provided geometry in the extended phase space", verbosity.low)

                if self.hessian_init =='true':
                    info(" @GEOP: We are computing the initial hessian", verbosity.low)
                    get_hessian(self.hessian, self.gm,self.beads.q)
                    add_spring_hessian(self.im, self.hessian)
                elif self.hessian_init =='expand':
                    # np.set_printoptions(precision=6, suppress=True)
                    expand_hessian(self.initial_hessian,self.hessian,self.beads.natoms)

                #Init im
                u,g = self.im(self.beads.q)


        self.old_x[:] = self.beads.q
        self.old_u[:] = self.forces.pot
        self.old_f[:] = self.forces.f

        if len(self.fixatoms) > 0:
            for dqb in self.old_f:
                dqb[self.fixatoms * 3] = 0.0
                dqb[self.fixatoms * 3 + 1] = 0.0
                dqb[self.fixatoms * 3 + 2] = 0.0


        # Do one step. Update hessian for the new position
        Instanton(self.old_x, self.old_f,self.im.f, self.hessian,self.hessian_update,self.hessian_asr, self.im,self.gm, self.big_step)

        # Update positions and forces
        self.beads.q = self.gm.dbeads.q
        self.forces.transfer_forces(self.gm.dforces)  # This forces the update of the forces

        # Exit simulation step
        d_x_max = np.amax(np.absolute(np.subtract(self.beads.q, self.old_x)))
        self.exitstep(self.forces.pot, self.old_u,d_x_max)

    def exitstep(self, fx, u0, x):

        #Modification of exitstep !ALBERTO spring energy
        """ Exits the simulation step. Computes time, checks for convergence. """
        info(" @GEOP: Updating bead positions", verbosity.debug)

        self.qtime += time.time()

        info(' @Exit step: Energy difference: %.1e, (condition: %.1e)' % (np.absolute((fx - u0) / self.beads.natoms)[0],self.tolerances["energy"] ),verbosity.medium)
        info(' @Exit step: Maximum force component: %.1e, (condition: %.1e)' % (np.amax(np.absolute(self.forces.f+self.im.f)), self.tolerances["force"]), verbosity.medium)
        info(' @Exit step: Maximum component step component: %.1e, (condition: %.1e)' % (x, self.tolerances["position"]), verbosity.medium)

        if (np.absolute((fx - u0) / self.beads.natoms) <= self.tolerances["energy"]) \
                and ((np.amax(np.absolute(self.forces.f + self.im.f)) <= self.tolerances["force"]) or
                         (np.linalg.norm(self.forces.f.flatten() - self.old_f.flatten()) <= 1e-08)) \
                and (x <= self.tolerances["position"]):

            print_instanton(self.hessian, self.gm,self.im,self.hessian_asr, self.final_rates)
            softexit.trigger("Geometry optimization converged. Exiting simulation")

        #Jeremy's condition
        #if (x <= self.tolerances["position"]):
        #    print_instanton(self.hessian, self.gm,self.im,self.hessian_asr, self.final_rates)
        #    softexit.trigger("Geometry optimization converged. Exiting simulation")

#-------------------------------------------------------------------------------------------------------------


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
    info(" @Instanton: Current step norm = %g" % d_x_norm, verbosity.medium)
    if d_x_norm > big_step:
        info(" @Instanton: Attempted step norm = %g, scaled down to %g" % (d_x_norm, big_step), verbosity.medium)
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
    elif update == 'recompute':
        get_hessian(h0,gm,x)
        if gm.dbeads.nbeads != 1:
            add_spring_hessian(im,h0)
        u, g = gm(x) #To update pos and for values in the mapper. (Yes, it is stricly unnecesary but keeps the code simple)


#-------------------------------------------------------------------------------------------------------------
def get_hessian(h,gm,x0,d=0.001):

#Think about the case you have numerical gradients

    info(" @Instanton: Computing hessian" ,verbosity.low)
    ii = gm.dbeads.natoms * 3
    h[:]=np.zeros((h.shape),float)

    for j in range(ii):
        info(" @Instanton: Computing hessian: %d of %d" % ((j+1),ii), verbosity.medium)
        x = x0.copy()

        x[:, j] = x0[:, j] + d
        e, f1 = gm(x)

        x[:, j] = x0[:, j] - d
        e, f2 = gm(x)
        g = (f1 - f2) / (2 * d)

        for i in range(gm.dbeads.nbeads):
            h[j + i * ii, i * ii:(i + 1) * ii] = g[i, :]

    gm.set_pos(x0)


def add_spring_hessian(im,h):
    """ Add spring terms to the extended hessian
        """
    if im.dbeads.nbeads == 1:
        return

    ii = im.dbeads.natoms * 3
    h_sp = im.dbeads.m3[0] * im.omega2

    # Diagonal
    diag1 = np.diag(h_sp)
    diag2 = np.diag(2.0 * h_sp)

    if im.mode == 'half':
        i=0
        h[i * ii:(i + 1) * ii, i * ii:(i + 1) * ii] += diag1
        i=im.dbeads.nbeads-1
        h[i * ii:(i + 1) * ii, i * ii:(i + 1) * ii] += diag1
        for i in range(1, im.dbeads.nbeads-1):
            h[i * ii:(i + 1) * ii, i * ii:(i + 1) * ii] += diag2
    elif im.mode =='full':
        for i in range(0, im.dbeads.nbeads):
            h[i * ii:(i + 1) * ii, i * ii:(i + 1) * ii] += diag2

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

    if h0 == None:
        return
    if h0.size != h.size:
        raise ValueError("@Expand hessian. This part is not ready yet.")

def get_imvector(h,  m3):
    """ Compute eigenvector  corresponding to the imaginary mode
            IN     h      = hessian
                   m3     = mass vector (dimension = 1 x 3*natoms)
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
        raise ValueError(" @GEOP: Small negative frequency %4.1f cm^-1" % freq, verbosity.low)
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
            q      = positions
            natoms = number of atoms
            nbeads = number of beads
            m      = mass vector, one value for each atom
            m3     = mass vector, one value for each degree of freedom
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
    #print dd[0:6]
    #Zeros
    condition = np.abs(dd) < 0.01 #Note that dd[] units are cm^1
    nzero= np.extract(condition, dd)

    if asr =='poly' and nzero.size != 6:
         info(" @GEOP: Warning, we have %d 'zero' frequencies" %nzero.size, verbosity.low)

    if asr =='crystal' and nzero.size != 3:
         info(" @GEOP: Warning, we have %d 'zero' frequencies" %nzero.size, verbosity.low)

    #Negatives
    condition = dd < -4.0 #Note that dd[] units are cm^1
    nneg = np.extract(condition, dd)
    info(" @Clean hessian: We have %d 'neg' frequencies " % (nneg.size), verbosity.medium)


    # Now eliminate external degrees of freedom from the dynmatrix
    print nzero.size
    if nzero.size > 0:
        if np.linalg.norm(nzero) > 0.01:
            info(" Warning @CLean hessian: We have deleted %d 'zero' frequencies " % (nzero.size), verbosity.high)
            info(" but the norm is greater than 0.01 cm^-1.  This should not happen." % (nzero.size), verbosity.high)

        d = np.delete(d,range(nneg.size,nneg.size+nzero.size))
        w = np.delete(w, range(nneg.size,nneg.size+nzero.size),axis=1)

    #dd = np.sign(d) * np.absolute(d) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17)
    #print dd[0:6]
    return d,w

def print_instanton(h,gm,im,asr,rates):
    # Compute action
    np.set_printoptions(precision=4, suppress=True,threshold=np.nan)
    action1 = gm.dforces.pot / (im.temp * im.dbeads.nbeads * units.Constants.kb)
    action2 = im.pot / (im.temp * im.dbeads.nbeads * units.Constants.kb)
    # Note that for the half polymer the factor 2 cancels out (One factor in *.pot and one factor in *.nbeads)
    print  'S1/hbar', action1 / units.Constants.hbar
    print  'S2/hbar', action2 / units.Constants.hbar
    print  'S/hbar', (action1 + action2) / units.Constants.hbar
    print  'ALBERTO', action1 / units.Constants.hbar, action2 / units.Constants.hbar, (action1 + action2) / units.Constants.hbar

    if rates != 'true':
        info(" Compute rates is false, we are not going to compute them.", verbosity.low)
    else:
        info(" Computing rates. For this we need to compute the hessian.", verbosity.low)


        if im.mode =='full':
            get_hessian(h, gm, im.dbeads.q)
            if gm.dbeads.nbeads != 1:
                add_spring_hessian(im, h)
            d,w = clean_hessian(h,im.dbeads.q,im.dbeads.natoms,im.dbeads.nbeads,im.dbeads.m,im.dbeads.m3,asr)
            print "Final  lowest eight frequencies"
            print np.sign(d[0:8]) * np.absolute(d[0:8]) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17)  # convert to cm^-1
        elif im.mode =='half':
            get_hessian(h, gm, im.dbeads.q)
            hbig=get_doble_hessian(h,im)
            q,nbeads,m,m3 = get_doble(im.dbeads.q, im.dbeads.nbeads, im.dbeads.m, im.dbeads.m3)
            d, w = clean_hessian(hbig, q, im.dbeads.natoms, nbeads, m, m3, asr)
            print "Final  lowest eight frequencies"
            print np.sign(d[0:8]) * np.absolute(d[0:8]) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17)  # convert to cm^-1

        compute_rates(im,gm)

def get_doble_hessian(h0,im):
    """Takes a hessian of the half polymer (only the physical part) and construct the 
       hessian for the full ringpolymer (physical + spring)"""

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

    ii = im.dbeads.natoms * 3

    # Spring part
    h_sp = im.dbeads.m3[0] * im.omega2

    # Diagonal
    diag = np.diag(2.0 * h_sp)

    for i in range(0, 2*im.dbeads.nbeads):
        h[i * ii:(i + 1) * ii, i * ii:(i + 1) * ii] += diag

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

def compute_rates(im,gm):
    info(" @GEOP: We are going to do the final postexi " , verbosity.low)
    pass
