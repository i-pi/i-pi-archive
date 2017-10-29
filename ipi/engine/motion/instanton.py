"""
Contains classes for instanton  calculations.

Algorithms implemented by Yair Litman and Mariana Rossi, 2017
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np
import time
import sys

from ipi.engine.motion import Motion
from ipi.utils.depend import depstrip, dobject
from ipi.utils.softexit import softexit
from ipi.utils.messages import verbosity, info
from ipi.utils import units
from ipi.utils.mintools import nichols, Powell
import scipy.sparse as sp
import scipy.sparse.linalg as sp_linalg

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
        ALBERTO
    """

    def __init__(self, fixcom=False, fixatoms=None,
                 mode="instanton",
                 tolerances={"energy": 1e-5, "force": 1e-5, "position": 5e-3},
                 biggest_step=0.3,
                 old_pos=np.zeros(0, float),
                 old_pot= np.zeros(0, float),
                 old_force=np.zeros(0, float),
                 opt='None',
                 action=np.zeros(2, float),
                 prefix="INSTANTON",
                 delta=np.zeros(0, float),
                 hessian_sparse=None,
                 hessian_init=None,
                 hessian=np.eye(0, 0, 0, float),
                 hessian_update=None,
                 hessian_asr=None,
                 qlist_lbfgs=np.zeros(0, float),
                 glist_lbfgs = np.zeros(0, float),
                 scale_lbfgs = 1,
                 corrections_lbfgs = 5,
                 final_post = 'False'):

        """Initialises InstantonMotion.
        """

        super(InstantonMotion, self).__init__(fixcom=fixcom, fixatoms=fixatoms)

        # Optimization mode
        self.mode           = mode

        # Generic optimization
        self.big_step       = biggest_step
        self.tolerances     = tolerances
        self.opt            = opt
        self.old_x          = old_pos
        self.old_u          = old_pot
        self.old_f          = old_force

        # Generic instanton
        self.action         = action
        self.prefix         = prefix
        self.delta          = delta
        self.final_post     = final_post

        #
        if mode =='rate':
            self.optimizer      = RateOptimizer()
            self.sparse         = hessian_sparse
            self.hessian_init   = hessian_init
            self.hessian        = hessian
            self.hessian_update = hessian_update
            self.hessian_asr = hessian_asr
        elif mode == 'splitting':
            # self.optimizer    = LBFGSOptimizer()
            self.corrections    = corrections_lbfgs
            self.scale          = scale_lbfgs
            self.qlist          = qlist_lbfgs
            self.glist          = glist_lbfgs


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

    """Creation of the multi-dimensional function to compute the physical potential and forces

    Attributes:
        dbeads:  copy of the bead object
        dcell:   copy of the cell object
        dforces: copy of the forces object
    """

    def __init__(self):
        pass

    def bind(self, dumop):
        self.dbeads  = dumop.beads.copy()
        self.dcell   = dumop.cell.copy()
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
    """Creation of the multi-dimensional function to compute full or half ring polymer pot
       and forces.
    """

    def __init__(self):
        self.pot = None
        self.f   = None
        pass

    def bind(self, dumop):
        self.dbeads = dumop.beads.copy()
        self.temp   = dumop.temp
        self.omega2 = (self.temp * (2*self.dbeads.nbeads) * units.Constants.kb / units.Constants.hbar) ** 2

        if dumop.sparse == 'false':
            self.h = self.spring_hessian()
        else:
            self.h_sparse = self.spring_sparse_hessian()

    def save(self,e,g):
        self.pot = e
        self.f   = -g

    def spring_hessian(self, mode='half'):
        """Compute the 'spring hessian'           
           IN     im      = instanton mapper

           OUT    h       = hessian with only the spring terms ('spring hessian')
            """

        #ALBERTO This is some how repited in the get_doble_hessian. Clean this.
        info(" @spring_hessian", verbosity.high)
        ii = self.dbeads.natoms * 3
        h = np.zeros([ii * self.dbeads.nbeads, ii * self.dbeads.nbeads])

        if self.dbeads.nbeads == 1:
            return h

        # Diagonal
        h_sp = self.dbeads.m3[0] * self.omega2
        diag1 = np.diag(h_sp)
        diag2 = np.diag(2.0 * h_sp)

        if mode == 'half':
            i = 0
            h[i * ii:(i + 1) * ii, i * ii:(i + 1) * ii] += diag1
            i = self.dbeads.nbeads - 1
            h[i * ii:(i + 1) * ii, i * ii:(i + 1) * ii] += diag1
            for i in range(1, self.dbeads.nbeads - 1):
                h[i * ii:(i + 1) * ii, i * ii:(i + 1) * ii] += diag2
        elif mode == 'full':
            for i in range(0, self.dbeads.nbeads):
                h[i * ii:(i + 1) * ii, i * ii:(i + 1) * ii] += diag2

        # Non-Diagonal
        ndiag = np.diag(-h_sp)
        # quasi-band
        for i in range(0, self.dbeads.nbeads - 1):
            h[i * ii:(i + 1) * ii, (i + 1) * ii:(i + 2) * ii] += ndiag
            h[(i + 1) * ii:(i + 2) * ii, i * ii:(i + 1) * ii] += ndiag

        # Corner
        if mode == 'full':
            h[0:ii, (self.dbeads.nbeads - 1) * ii:(self.dbeads.nbeads) * ii] += ndiag
            h[(self.dbeads.nbeads - 1) * ii:(self.dbeads.nbeads) * ii, 0:ii] += ndiag

        return h

    def spring_sparse_hessian(self):
        """Create the spring hessian sparse matrix for the half RP'           
           IN     im      = instanton mapper

           OUT    h       = sparse hessian with only the spring terms ('spring hessian')
            """

        info(" @spring_hessian_sparse", verbosity.high)

        ii = self.dbeads.natoms * 3


        # Diagonal
        d_corner  = self.dbeads.m3[0] * self.omega2
        d_0       = np.array([[d_corner*2]]).repeat(self.dbeads.nbeads-2,axis=0).flatten()
        diag_0    = np.concatenate((d_corner,d_0,d_corner))

        # Non-Diagonal
        d_out     = - d_corner
        diag_ii   = np.array([[d_out]]).repeat(self.dbeads.nbeads,axis=0).flatten()

        #Create matrix
        offsets   = np.array([0, -ii, ii])
        h         = sp.dia_matrix((np.array([diag_0,diag_ii,diag_ii]),offsets),shape=(self.dbeads.q.size,self.dbeads.q.size))

        return h



    def __call__(self, x,ret=True):
        """Computes spring energy and gradient for instanton optimization step"""

        if x.shape[0]==1: #only one bead
           self.dbeads.q = x.copy()
           e=0.0
           g=np.zeros(x.shape[1])
           self.save(e, g)
        else:
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

        if ret:
            return e, g


class DummyOptimizer(dobject):
    """ Dummy class for all optimization classes """

    def __init__(self):

        """Initialises object for GradientMapper (physical potential, forces and hessian) 
        and InstantonMapper ( spring potential,forces and hessian) """

        self.gm           = GradientMapper()
        self.im           = InstantonMapper()
        self.exit         = False

    def step(self, step=None):
        """Dummy simulation time step which does nothing."""
        pass

    def bind(self, geop):

        """
        bind optimization options and call bind function of Mappers (get beads, cell,forces)
        check whether force size,  Hessian size from  match system size
        """
        if geop.beads.nbeads == 1 and geop.mode=='half':
            raise ValueError("Half mode is not compatible with nbeads = 1")


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

        #Temperature
        self.temp = geop.ensemble.temp
        if geop.ensemble.temp == -1.0 or geop.ensemble.temp == 1.0:  # This is due to a little inconsistency on the default value
            if self.beads.nbeads != 1:
                raise ValueError("Temperature must be specified for an Instanton calculation ")

        # Optimization mode
        self.mode       = geop.mode

        # Generic optimization
        self.tolerances = geop.tolerances
        self.big_step   = geop.big_step
        self.old_x      = geop.old_x
        self.old_u      = geop.old_u
        self.old_f      = geop.old_f

        # Generic instanton

        if geop.action.size !=2:
            geop.action      = np.zeros(2, float)
        self.action          = geop.action
        self.prefix          = geop.prefix
        self.delta           = geop.delta
        self.final_post      = geop.final_post
        self.gm.bind(self)

    def exitstep(self, fx, fx0, x,exitt):

        """ Exits the simulation step. Computes time, checks for convergence. """
        self.qtime += time.time()

        info(' @Exit step: Energy difference: %.1e, (condition: %.1e)' % (np.absolute((fx - fx0) / self.beads.natoms)[0],self.tolerances["energy"] ),verbosity.low)
        info(' @Exit step: Maximum force component: %.1e, (condition: %.1e)' % (np.amax(np.absolute(self.forces.f+self.im.f)), self.tolerances["force"]), verbosity.low)
        info(' @Exit step: Maximum component step component: %.1e, (condition: %.1e)' % (x, self.tolerances["position"]), verbosity.low)

        if (np.absolute((fx - fx0) / self.beads.natoms) <= self.tolerances["energy"]) \
                and ((np.amax(np.absolute(self.forces.f + self.im.f)) <= self.tolerances["force"]) or
                         (np.linalg.norm(self.forces.f.flatten() - self.old_f.flatten()) <= 1e-08)) \
                and (x <= self.tolerances["position"]):

            print_instanton(self.prefix, self.hessian, self.gm,self.im,self.hessian_asr,self.final_post,self.action)
            exitt=True #If we just exit here, the last step (including the last hessian) will not be in the RESTART file

        return exitt

class RateOptimizer(DummyOptimizer):
    """ INSTANTON Rate calculation"""

    def bind(self, geop):
        # call bind function from DummyOptimizer
        super(RateOptimizer,self).bind(geop)

        #Specific for RateOptimizer
        self.sparse          = geop.sparse
        self.hessian_update  = geop.hessian_update
        self.hessian_asr     = geop.hessian_asr
        self.hessian_init    = geop.hessian_init

        self.im.bind(self)

        #Hessian
        self.initial_hessian = None

        if geop.hessian.size != (self.beads.natoms*3 * self.beads.q.size ):
            if geop.hessian.size == (self.beads.natoms*3)**2:
                self.initial_hessian = geop.hessian.copy()
                geop.hessian = np.zeros((self.beads.natoms*3, self.beads.q.size), float)
            elif self.beads.nbeads == 1:
               if geop.hessian.size == 0 and geop.hessian_init == 'true':
                   info(" Initial classical hessian is not provided. We are going to compute it.", verbosity.low)
                   geop.hessian = np.zeros((self.beads.natoms*3, self.beads.q.size))
               else:
                   raise ValueError("Nbeads =1. Hessian_init != 'True'. An initial hessian (size natoms*3 X natoms*3) must be provided")
            else:
                raise ValueError("Nbeads >1. An initial hessian (shape '(natoms*3,natoms*3)' or '(natoms*3,natoms*nbeads*3)') must be provided.\n Hessian size does not match system size.")

        self.hessian = geop.hessian


    def step(self, step=None):
        """ Does one simulation time step."""

        self.qtime = -time.time()
        info("\n Instanton optimization STEP %d" % step, verbosity.low)


        if step == 0:
            info(" @GEOP: Initializing INSTANTON", verbosity.low)

            if self.beads.nbeads == 1:
                info(" @GEOP: Classical TS search", verbosity.low)
                if self.hessian_init == 'true':
                    get_hessian(self.hessian, self.gm, self.beads.q)
            else:
                if ((self.beads.q - self.beads.q[0]) == 0).all():  # If the coordinates in all the imaginary time slices are the same
                    info(" @GEOP: We stretch the initial geometry with an 'amplitud' of %4.2f" % self.delta, verbosity.low)
                    imvector = get_imvector(self.initial_hessian,self.beads.m3[0].flatten())
                    for i in range(self.beads.nbeads):
                        self.beads.q[i, :] += self.delta * np.cos(i * np.pi / float(self.beads.nbeads-1)) * imvector[:]
                    if self.hessian_init !='true':
                        info(" @GEOP: Hessian_init isn't true but we have stretched the polymer so we are going to compute the initial hessian anyway.", verbosity.low)
                    self.hessian_init ='true'
                else:
                    info(" @GEOP: Starting from the provided geometry in the extended phase space", verbosity.low)

                if self.hessian_init =='true':
                    info(" @GEOP: We are computing the initial hessian", verbosity.low)
                    get_hessian(self.hessian, self.gm,self.beads.q)

        if self.exit:
            softexit.trigger("Geometry optimization converged. Exiting simulation")

        if self.im.f== None:
            u,g = self.im(self.beads.q) #Init instanton mapper


        self.old_x[:] = self.beads.q
        self.old_u[:] = self.forces.pot
        self.old_f[:] = self.forces.f

        if len(self.fixatoms) > 0:
            for dqb in self.old_f:
                dqb[self.fixatoms * 3] = 0.0
                dqb[self.fixatoms * 3 + 1] = 0.0
                dqb[self.fixatoms * 3 + 2] = 0.0


        # Do one step. Update hessian for the new position. Update the position and force inside the mapper.
        Instanton(self.old_x, self.old_f,self.im.f, self.sparse,self.hessian,self.hessian_update,self.action,self.hessian_asr, self.im,self.gm, self.big_step)

        # Update positions and forces
        self.beads.q = self.gm.dbeads.q
        self.forces.transfer_forces(self.gm.dforces)  # This forces the update of the forces

        # Exit simulation step
        d_x_max = np.amax(np.absolute(np.subtract(self.beads.q, self.old_x)))
        self.exit=self.exitstep(self.forces.pot, self.old_u,d_x_max,self.exit)



#-------------------------------------------------------------------------------------------------------------
def Instanton(x0, f0,f1, sparse,h, update,action,asr, im,gm, big_step):
    """Do one step. Update hessian for the new position. Update the position and force inside the mapper.
       
       Input: x0  = last positions
               f0 = last physical forces
               f1 = last spring forces
               h  = physical hessian
           update = how to update the hessian
           action = vector to store the current value of the action
               im = instanton mapper
               gm = gradient  mapper
         big_step = limit on step length"""

    info(" @Instanton_step", verbosity.high)

    #Construct hessian and get eigenvalues and eigenvector

    if sparse == 'false':
        h0      = red2comp(h,im)    #construct complete hessian from reduced
        h1      = np.add(im.h, h0)  # add spring terms to the physical hessian
        d, w    = clean_hessian(h1, im.dbeads.q, im.dbeads.natoms, im.dbeads.nbeads, im.dbeads.m, im.dbeads.m3, asr)
    else:
        h0      = red2sparse(h)        #construct sparse hessian from reduced
        h1      = ( h0 + im.h_sparse ) # add spring terms to the physical hessian
        d, w    = clean_sparse_hessian(h1, im.dbeads.q, im.dbeads.natoms, im.dbeads.nbeads, im.dbeads.m, im.dbeads.m3,asr)


    # Find new movement direction
    d_x = nichols(f0,f1, d,w,im.dbeads.m3, big_step)

    # Rescale step
    d_x_max = np.amax(np.absolute(d_x))
    info(" @Instanton: Current step norm = %g" % d_x_max, verbosity.medium)
    if np.amax(np.absolute(d_x)) > big_step:
        info(" @Instanton: Attempted step norm = %g, scaled down to %g" % (d_x_max, big_step), verbosity.low)
        d_x *= big_step / np.amax(np.absolute(d_x_max))

    # Make movement and get new energy (u)  and forces(f) using mapper
    x = x0 + d_x
    im(x,ret=False) # Only to update the mapper
    u, g2 = gm(x)
    f = -g2

    # Update hessian
    if update == 'powell':
        d_g = np.subtract(f0, f)
        ##Powell(d_x.flatten(), d_g.flatten(), h0)
        ##comp2red(h, h0, im)
        i = im.dbeads.natoms * 3
        for j in range(im.dbeads.nbeads):
            aux = h[:, j * i:(j + 1) * i]
            dg = d_g[j, :]
            dx = d_x[j, :]
            Powell(dx, dg, aux)

    elif update == 'recompute':
        get_hessian(h,gm,x)

    #Store action
    action[0] = gm.dforces.pot * 1/(im.temp  * im.dbeads.nbeads * units.Constants.kb)
    action[1] = im.pot / (im.temp * im.dbeads.nbeads * units.Constants.kb)
    # Note that for the half polymer the factor 2 cancels out (One factor in *.pot and one factor in *.nbeads)

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------

#Alberto. All this hessian related functions can be join inside a hessian object.
class hessian(object):
    def __init__(self):
        pass

#--------------------------------------------------------------------------------------------------------------
def red2comp(h, im):
    """Takes the reduced physical hessian and construct the 'complete' one (all 0 included) """
    info("\n @Instanton: Creating 'complete' physical hessian \n", verbosity.medium)
    i = im.dbeads.natoms * 3
    ii = im.dbeads.q.size
    h0 = np.zeros((ii, ii), float)

    for j in range(im.dbeads.nbeads):
        h0[j * i:(j + 1) * i, j * i:(j + 1) * i] = h[:, j * i:(j + 1) * i]
    return h0


def comp2red(h, h0, im):
    """Takes the 'complete' physical hessian and construct the reduced one  """

    i = im.dbeads.natoms * 3
    h[:] = np.zeros((h.shape), float)

    for j in range(im.dbeads.nbeads):
        h[:, j * i:(j + 1) * i] = h0[j * i:(j + 1) * i, j * i:(j + 1) * i]

#--------------------------------------------------------------------------------------------------------------
def red2sparse(h):
    i = h.shape[0]   #natoms*3
    n = h.shape[1]/i #nbeads
    aux = list()
    for j in range(n):
        aux.append(h[:,j*i:(j+1)*i])
    hsparse=sp.block_diag(aux,format='lil')
    #np.set_printoptions(precision=4, suppress=True,threshold=np.nan,linewidth=500)

    return hsparse

def get_hessian(h,gm,x0,d=0.01):
    """Compute the physical hessian           
       IN     h       = physical hessian 
              gm      = gradient mapper
              x0      = position vector 
              
       OUT    h       = physical hessian
        """
#TODO What about the case you have numerical gradients?


    info(" @Instanton: Computing hessian" ,verbosity.low)
    ii = gm.dbeads.natoms * 3
    h[:]=np.zeros((h.shape),float)

    ##Ask Michelle about transfer force here I
    #ddbeads = gm.dbeads.copy()
    #ddcell = gm.dcell.copy()
    #ddforces = gm.dforces.copy(ddbeads, ddcell)


    for j in range(ii):
        info(" @Instanton: Computing hessian: %d of %d" % ((j+1),ii), verbosity.low)
        x = x0.copy()

        x[:, j] = x0[:, j] + d
        e, f1 = gm(x)

        x[:, j] = x0[:, j] - d
        e, f2 = gm(x)
        g = (f1 - f2) / (2 * d)

        for i in range(gm.dbeads.nbeads):
            h[j, :] = g.flatten()

        #for i in range(gm.dbeads.nbeads):
        #    h[j + i * ii, i * ii:(i + 1) * ii] = g[i, :]

    u,g=gm(x0) #Keep the mapper updated

    ## Ask Michelle about transfer force here II
    #gm.dbeads.q = ddbeads.q
    #gm.dforces.transfer_forces(ddforces)


def clean_hessian(h, q, natoms, nbeads, m, m3, asr, mofi=False):
    """
        Removes the translations and rotations modes.
        IN  h      = hessian
            q      = positions
            natoms = number of atoms
            nbeads = number of beads
            m      = mass vector, one value for each atom
            m3     = mass vector, one value for each degree of freedom
            mofi   = An optional boolean which decides whether the det(M_of_I)
                     is returned or not. Defaults to False.
        OUT d      = non zero eigenvalues of the dynmatrix
            w      = eigenvectors without external modes

        #Adapted from ipi/engine/motion/phonons.py apply_asr    """

    info(" @clean_hessian", verbosity.high)
    # Set some useful things
    ii = natoms * nbeads
    mm = np.zeros((nbeads, natoms))
    for i in range(nbeads):
        mm[i] = m
    mm = mm.reshape(ii)
    ism = m3.reshape(ii * 3) ** (-0.5)
    ismm = np.outer(ism, ism)
    dynmat = np.multiply(h, ismm)

    if asr == 'none':
        hm = dynmat
    else:
        # Computes the centre of mass.
        com = np.dot(np.transpose(q.reshape((ii, 3))), mm) / mm.sum()
        qminuscom = q.reshape((ii, 3)) - com

        if asr == 'poly':
            # Computes the moment of inertia tensor.
            moi = np.zeros((3, 3), float)
            for k in range(ii):
                moi -= np.dot(np.cross(qminuscom[k], np.identity(3)), np.cross(qminuscom[k], np.identity(3))) * mm[k]

            I, U = (np.linalg.eig(moi))
            R = np.dot(qminuscom, U)
            D = np.zeros((6, 3 * ii), float)

            # Computes the vectors along translations and rotations.
            # Translations
            D[0] = np.tile([1, 0, 0], ii) / ism
            D[1] = np.tile([0, 1, 0], ii) / ism
            D[2] = np.tile([0, 0, 1], ii) / ism
            # Rotations
            for i in range(3 * ii):
                iatom = i / 3
                idof = np.mod(i, 3)
                D[3, i] = (R[iatom, 1] * U[idof, 2] - R[iatom, 2] * U[idof, 1]) / ism[i]
                D[4, i] = (R[iatom, 2] * U[idof, 0] - R[iatom, 0] * U[idof, 2]) / ism[i]
                D[5, i] = (R[iatom, 0] * U[idof, 1] - R[iatom, 1] * U[idof, 0]) / ism[i]

            for k in range(6):
                D[k] = D[k] / np.linalg.norm(D[k])
            # Computes the transformation matrix.
            transfmatrix = np.eye(3 * ii) - np.dot(D.T, D)
            hm = np.dot(transfmatrix.T, np.dot(dynmat, transfmatrix))

        elif asr == 'crystal':
            # Computes the vectors along translations.
            # Translations
            D = np.zeros((3, 3 * ii), float)
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
    hm = (hmT + hm) / 2.0

    d, w = np.linalg.eigh(hm)

    # Count
    dd = np.sign(d) * np.absolute(d) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17)  # convert to cm^-1
    #print dd[0:9]
    # Zeros
    cut0 = 0.01  # Note that dd[] units are cm^1
    condition = np.abs(dd) < cut0
    nzero = np.extract(condition, dd)

    if asr == 'poly' and nzero.size != 6:
        info(" @GEOP: Warning, we have %d 'zero' frequencies" % nzero.size, verbosity.low)

    if asr == 'crystal' and nzero.size != 3:
        info(" @GEOP: Warning, we have %d 'zero' frequencies" % nzero.size, verbosity.low)

    # Negatives
    cutNeg = -4  # Note that dd[] units are cm^1
    condition = dd < cutNeg
    nneg = np.extract(condition, dd)
    info(" @Clean hessian: We have %d 'neg' frequencies " % (nneg.size), verbosity.medium)

    # Now eliminate external degrees of freedom from the dynmatrix

    if nzero.size > 0:
        if np.linalg.norm(nzero) > cut0:
            info(" Warning @Clean hessian: We have deleted %d 'zero' frequencies " % (nzero.size), verbosity.high)
            info(" but the norm is greater than 0.01 cm^-1.  This should not happen.", verbosity.high)

        d = np.delete(d, range(nneg.size, nneg.size + nzero.size))
        w = np.delete(w, range(nneg.size, nneg.size + nzero.size), axis=1)

    if mofi:
        if asr == 'poly':
            return d, w, np.prod(I)
        else:
            return d, w, 1.0
    else:
        return d, w

# ----------------------------------------------------------------------------------------------------------

def clean_sparse_hessian(h, q, natoms, nbeads, m, m3, asr, mofi=False):
    """
        Removes the translations and rotations modes.
        IN  h      = sparse hessian
            q      = positions
            natoms = number of atoms
            nbeads = number of beads
            m      = mass vector, one value for each atom
            m3     = mass vector, one value for each degree of freedom
            mofi   = An optional boolean which decides whether the det(M_of_I)
                     is returned or not. Defaults to False.
        OUT d      = non zero eigenvalues of the dynmatrix
            w      = dynmatrix with the external modes projected out

        #Adapted from ipi/engine/motion/phonons.py apply_asr    """
    info(" @clean_sparse_hessian", verbosity.high)

    # Set some useful things
    ii = natoms * nbeads
    ism = m3.reshape(ii * 3) ** (-0.5)
    diag = np.array([ism]).repeat(nbeads, axis=0).flatten()
    offsets = np.array([0])

    ISM = sp.dia_matrix((np.array([diag]), offsets), shape=(q.size,q.size))

    dynmat=ISM.dot(h.dot(ISM))

    ##Simmetrize to use linalg.eigh
    hm  = dynmat

    hmT = hm.transpose(copy=True)
    hmm = (hmT + hm) / 2.0

    #Work HERE
    d1, w1 = sp_linalg.eigsh(hmm, which='SA', k=ii*3/2)  #Smallest k
    d2, w2 = sp_linalg.eigsh(hmm, which='LA', k=ii*3/2) #Largest k
    d = np.concatenate((d1,d2))
    w = np.concatenate((w1, w2),axis=1)

    # Count
    dd = np.sign(d) * np.absolute(d) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17)  # convert to cm^-1
    np.set_printoptions(precision=6, suppress=True, threshold=np.nan)
    print dd[0:9]

    # Zeros
    cut0 = 30  # Note that dd[] units are cm^1
    condition = np.abs(dd) < cut0
    nzero = np.extract(condition, dd)

    # Negatives
    cutNeg = -30  # Note that dd[] units are cm^1
    condition = dd < cutNeg
    nneg = np.extract(condition, dd)

    info(" @Clean hessian: We have %d 'neg' frequencies " % (nneg.size), verbosity.medium)

    # Now eliminate external degrees of freedom from the dynmatrix

    if nzero.size > 0:
        if np.linalg.norm(nzero) > cut0:
            info(" Warning @Clean hessian: We have deleted %d 'zero' frequencies " % (nzero.size), verbosity.high)
            info(" but the norm is greater than 0.01 cm^-1.  This should not happen.", verbosity.high)

        d = np.delete(d, range(nneg.size, nneg.size + nzero.size))
        w = np.delete(w, range(nneg.size, nneg.size + nzero.size), axis=1)

    if mofi:
        if asr == 'poly':
            #This is repited. Clean this.
            # Computes the moment of inertia tensor.
            mm = np.zeros((nbeads, natoms))
            for i in range(nbeads):
                mm[i] = m
            mm = mm.reshape(ii)
            com = np.dot(np.transpose(q.reshape((ii, 3))), mm) / mm.sum()
            qminuscom = q.reshape((ii, 3)) - com
            moi = np.zeros((3, 3), float)
            for k in range(ii):
                moi -= np.dot(np.cross(qminuscom[k], np.identity(3)), np.cross(qminuscom[k], np.identity(3))) * mm[k]

            I, U = (np.linalg.eig(moi))

            return d, w, np.prod(I)
        else:
            return d, w, 1.0
    else:
        return d, w



# ----------------------------------------------------------------------------------------------------------

def get_doble_hessian(h0, im):
    """Takes a hessian of the half polymer (only the physical part) and construct the 
       hessian for the full polymer (physical + spring).

       IN  h0      = physical hessian of the half polymer
            im     = instanton mapper

       OUT h      = full ring polymer hessian (physical + spring terms)"""
    info(" @get_doble_hessian", verbosity.high)

    nbeads = im.dbeads.nbeads
    natoms = im.dbeads.natoms
    ii = 3 * natoms
    iii = 3 * natoms * nbeads

    # np.set_printoptions(precision=4, suppress=True)

    h = np.zeros((iii * 2, iii * 2))
    h[0:iii, 0:iii] = h0

    # diagonal block
    for i in range(nbeads):
        x = i * ii + iii
        y = ((nbeads - 1) - i) * ii
        h[x:x + ii, x:x + ii] = h0[y:y + ii, y:y + ii]

    # Spring part
    h_sp = im.dbeads.m3[0] * im.omega2

    # Diagonal
    diag = np.diag(2.0 * h_sp)

    for i in range(0, 2 * im.dbeads.nbeads):
        h[i * ii:(i + 1) * ii, i * ii:(i + 1) * ii] += diag

    # quasi-band
    ndiag = np.diag(-h_sp)
    for i in range(0, 2 * im.dbeads.nbeads - 1):
        h[i * ii:(i + 1) * ii, (i + 1) * ii:(i + 2) * ii] += ndiag
        h[(i + 1) * ii:(i + 2) * ii, i * ii:(i + 1) * ii] += ndiag

    # Corner
    h[0:ii, (2 * nbeads - 1) * ii:(2 * nbeads) * ii] += ndiag
    h[(2 * nbeads - 1) * ii:(2 * nbeads) * ii, 0:ii] += ndiag

    return h

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------


def get_imvector(h,  m3):
    """ Compute eigenvector  corresponding to the imaginary mode
            IN     h      = hessian
                   m3     = mass vector (dimension = 1 x 3*natoms)
            OUT    imv    = eigenvector corresponding to the imaginary mode
        """
    info("@get_imvector", verbosity.high)
    if h.size != m3.size**2:
        raise ValueError("@Get_imvector. Initial hessian size does not match system size.")
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

    info(" @get_imvector: We stretch along the mode with freq %f cm^1" %freq[0] , verbosity.low)

    imv = w[:, 0] * (m3[:] ** 0.5)
    imv = imv/np.linalg.norm(imv)

    return imv

#-------------------------------------------------------------------------------------------------------------------

def get_doble(q0,nbeads,m,m3):
    """Takes a positions and mass vectors of the half polymer and construct the 
       equivalent for the full ringpolymer
       IN  q0     = positions of the half polymer
           nbeads = number of beads
           m      = mass vector
           m3     = extended mass vector, one entry per each degree of freedom
       OUT
           The corresponding vectors/values for the full ring-polymer
       """
    info(" @get_doble", verbosity.high)
    q=np.concatenate((q0, np.flipud(q0)), axis=0)
    m3=np.concatenate((m3, m3), axis=0)
    return q,2*nbeads,m,m3

#----------------------------------------------------------------------------------------------------------
def print_instanton(prefix,h,gm,im,asr,rates,action):


    """ Prints out relevant information."""
    info(" @print_instanton", verbosity.high)
    outfile = open(prefix + '.data', 'w')
    np.set_printoptions(precision=6, suppress=True, threshold=np.nan)

    #factor=2.0

    # action1 = gm.dforces.pot * 1/(im.temp  * im.dbeads.nbeads * units.Constants.kb)
    # action2 = im.pot / (im.temp * im.dbeads.nbeads * units.Constants.kb)
    # Note that for the half polymer the factor 2 cancels out (One factor in *.pot and one factor in *.nbeads)
    action1 = action[0]
    action2 = action[1]
    action = action1 + action2
    BN = 2 * np.sum(gm.dbeads.m3[1:, :] * (gm.dbeads.q[1:, :] - gm.dbeads.q[:-1, :]) ** 2)

    print >> outfile, "# Instanton calculation:"
    print >> outfile, ('  ')

    if im.dbeads.nbeads>1:
        print >> outfile, ('We have %i total beads (i.e. %i different) and the temperature  is %f K ' %(2*im.dbeads.nbeads,im.dbeads.nbeads,units.unit_to_user('temperature',"kelvin",im.temp)))
        print >> outfile, ('Beta in a.u.', 1 / im.temp)
        print >> outfile, ('S1/hbar   '  + str(action1 ) )
        print >> outfile, ('S2/hbar   '  + str(action2 ) )
        print >> outfile, ('S/hbar    '  + str(action  ) )
        print >> outfile, ('BN   ' + (str(BN)))
        print >> outfile, ('BN*N   ' + (str(BN * 2.0 * im.dbeads.nbeads)))
        print >> outfile, ('  ')
    else:
        print >> outfile, ('We have %i total beads. Classical TS search.' %im.dbeads.nbeads)
        print >> outfile, ('Potential (eV)   ' + str(units.unit_to_user('energy',"electronvolt",action1)*(im.temp  * im.dbeads.nbeads * units.Constants.kb)*1))
    info(" ", verbosity.low)
    info(" We are almost done.", verbosity.low)
    if rates != 'true':
        info(" Compute rates is false, we are not going to compute them.", verbosity.low)
    else:
        info(" Computing rates. For this we need to compute the hessian.", verbosity.low)

        get_hessian(h, gm, im.dbeads.q)

        if im.dbeads.nbeads>1:
            h0 = red2comp(h, im)
            hbig = get_doble_hessian(h0,im)
            q,nbeads,m,m3 = get_doble(im.dbeads.q, im.dbeads.nbeads, im.dbeads.m, im.dbeads.m3)
            d, w,detI = clean_hessian(hbig, q, im.dbeads.natoms, nbeads, m, m3, asr,mofi=True)
            print >> outfile, "Final  lowest ten frequencies (cm^-1)"
            print >> outfile, str(np.sign(d[0:10]) * np.absolute(d[0:10]) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17))  # convert to cm^-1
            print >> outfile, ('  ')

            # Compute Qrot
            # Compute Qtras

            # Compute loqQvib
            # if im.dbeads.nbeads ==1: # TS calculation,exact expression
            #    logQvib = np.sum(np.log(2 * np.sinh((betaP * units.Constants.hbar * np.sqrt(np.delete(d, 0)) / 2.0))))
            #    # logQvib = np.sum( np.log( betaP * units.Constants.hbar * np.sqrt(np.absolute(np.delete(d, 0)))))
            #    print 'Deleted frequency for computing Qvib  %f cm^-1' % (np.sign(d[0]) * np.absolute(d[0]) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17))
            # else: # instanton calculation
            #    logQvib = np.sum(np.log(betaP * units.Constants.hbar * np.sqrt(np.absolute(np.delete(d, 1)))))
            #    print 'Deleted frequency for computing Qvib  %f cm^-1' % ( np.sign(d[1]) * np.absolute(d[1]) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17))

            # print >> outfile, ('logQvib ' + str(logQvib))
        else: #TS
            d, w = clean_hessian(h, im.dbeads.q, im.dbeads.natoms, im.dbeads.nbeads, im.dbeads.m, im.dbeads.m3, asr)
            print >> outfile, "Frequencies (cm^-1)"
            print >> outfile, str( np.sign(d) * np.absolute(d) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17) ) # convert to cm^-1
            print >> outfile, ('  ')


    outfile.close()
