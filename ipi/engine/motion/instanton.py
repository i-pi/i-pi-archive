"""
Contains classes for instanton rates calculation.

Algorithms implemented by Yair Litman and Mariana Rossi, 2017
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np
import time

from ipi.engine.motion import Motion
from ipi.utils.depend import depstrip, dobject
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
        ALBERTO
    """

    def __init__(self, fixcom=False, fixatoms=None,
                 mode="instanton",
                 biggest_step=0.5,
                 old_pos=np.zeros(0, float),
                 old_pot= np.zeros(0, float),
                 old_force=np.zeros(0, float),
                 hessian=np.eye(0, 0, 0, float),
                 tolerances={"energy": 1e-5, "force": 1e-5, "position": 5e-3},
                 delta=np.zeros(0, float),
                 hessian_init=None,
                 hessian_update=None,
                 hessian_asr=None,
                 action=np.zeros(2, float),
                 prefix="INSTANTON",
                 final_rates='False'):

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
        self.prefix         = prefix
        self.action         = action
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

    """Creation of the multi-dimensional function to compute the physical potential and forces

    Attributes:
        dbeads:  copy of the bead object
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
        self.mode   = dumop.mode
        if self.mode=='full':
            self.omega2 = (self.temp * self.dbeads.nbeads * units.Constants.kb / units.Constants.hbar) ** 2
        elif self.mode =='half':
            self.omega2 = (self.temp * (2*self.dbeads.nbeads) * units.Constants.kb / units.Constants.hbar) ** 2
        self.h = spring_hessian(self)

    def save(self,e,g):
        self.pot = e
        self.f   = -g

    def __call__(self, x,ret=True):
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

        if ret:
            return e, g

class InstantonOptimizer(dobject):
    """ INSTANTON """

    def __init__(self):

        """Initialises object for GradientMapper (physical potential, forces and hessian) and InstantonMapper ( spring potential,forces and hessian) """

        self.gm           = GradientMapper()
        self.im           = InstantonMapper()
        self.exit         = False

    def bind(self, geop):

        """
        bind optimization options and call bind function of Mappers (get beads, cell,forces)
        check whether force size,  Hessian size from  match system size
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

        if geop.action.size !=2:
            geop.action = np.zeros(2, float)
        self.action          = geop.action
        self.prefix          = geop.prefix

        self.gm.bind(self)
        self.im.bind(self)

        #Hessian
        self.initial_hessian = None

        if geop.hessian.size != (self.beads.q.size * self.beads.q.size):
            if geop.hessian.size == (self.beads.natoms*3)**2:
                self.initial_hessian = geop.hessian.copy()
                geop.hessian = np.zeros((self.beads.q.size, self.beads.q.size), float)
            elif self.beads.nbeads == 1:
               if geop.hessian.size == 0 and geop.hessian_init == 'true':
                   info(" Initial classical hessian is not provided. We are going to compute it.", verbosity.low)
                   geop.hessian = np.zeros((self.beads.q.size, self.beads.q.size))
               else:
                   raise ValueError("Nbeads =1. Hessian_init != 'True'. An initial hessian (size natoms*3 X natoms*3) must be provided")
            else:
                raise ValueError("Nbeads >1. An initial hessian (size 'natoms*3^2' or 'natoms*nbeads*3^2') must be providedHessian size does not match system size.")

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
                        if self.mode == 'full':
                            self.beads.q[i, :] += self.delta * np.cos((i+1) *2.0*np.pi / float(self.beads.nbeads)) * imvector[:]
                        elif self.mode== 'half':
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
        Instanton(self.old_x, self.old_f,self.im.f, self.hessian,self.hessian_update,self.action,self.hessian_asr, self.im,self.gm, self.big_step)

        # Update positions and forces
        self.beads.q = self.gm.dbeads.q
        self.forces.transfer_forces(self.gm.dforces)  # This forces the update of the forces

        # Exit simulation step
        d_x_max = np.amax(np.absolute(np.subtract(self.beads.q, self.old_x)))
        self.exit=self.exitstep(self.forces.pot, self.old_u,d_x_max,self.exit)

    def exitstep(self, fx, fx0, x,exitt):

        #Modification of exitstep !ALBERTO spring energy
        """ Exits the simulation step. Computes time, checks for convergence. """
        self.qtime += time.time()

        info(' @Exit step: Energy difference: %.1e, (condition: %.1e)' % (np.absolute((fx - fx0) / self.beads.natoms)[0],self.tolerances["energy"] ),verbosity.low)
        info(' @Exit step: Maximum force component: %.1e, (condition: %.1e)' % (np.amax(np.absolute(self.forces.f+self.im.f)), self.tolerances["force"]), verbosity.low)
        info(' @Exit step: Maximum component step component: %.1e, (condition: %.1e)' % (x, self.tolerances["position"]), verbosity.low)

        if (np.absolute((fx - fx0) / self.beads.natoms) <= self.tolerances["energy"]) \
                and ((np.amax(np.absolute(self.forces.f + self.im.f)) <= self.tolerances["force"]) or
                         (np.linalg.norm(self.forces.f.flatten() - self.old_f.flatten()) <= 1e-08)) \
                and (x <= self.tolerances["position"]):

            print_instanton(self.prefix, self.hessian, self.gm,self.im,self.hessian_asr,self.final_rates,self.action)
            exitt=True #If we just exit here, the last step (including the last hessian) will not be in the RESTART file

        return exitt

#-------------------------------------------------------------------------------------------------------------


def Instanton(x0, f0,f1, h, update,action,asr, im,gm, big_step):
    """Do one step. Update hessian for the new position. Update the position and force inside the mapper.
       
       Input: x0 = last positions
               f0 = last physical forces
               f1 = last spring forces
               h  = physical hessian
           update = how to update the hessian
           action = vector to store the current value of the action
               im = instanton mapper
               gm = gradient  mapper
         big_step = limit on step length"""

    info(" @Instanton_step", verbosity.high)

    # Project out rotations and translation from the Hessian. Note that the dynmax.size < h0
    time0 = time.time()
    h1 = np.add(im.h, h) #add spring terms to the physical hessian
    time1 = time.time()
    d, dynmax = clean_hessian(h1, im.dbeads.q, im.dbeads.natoms, im.dbeads.nbeads, im.dbeads.m, im.dbeads.m3,asr)


    # Find new movement direction
    time2 = time.time()
    d_x = nichols(f0,f1, d,dynmax,im.dbeads.m3, big_step)

    # Rescale step
    time3 = time.time()
    d_x_max = np.amax(np.absolute(d_x))
    info(" @Instanton: Current step norm = %g" % d_x_max, verbosity.medium)
    if np.amax(np.absolute(d_x)) > big_step:
        info(" @Instanton: Attempted step norm = %g, scaled down to %g" % (d_x_max, big_step), verbosity.low)
        d_x *= big_step / np.amax(np.absolute(d_x_max))

    # Make movement and get new energy (u)  and forces(f) using mapper
    time4 = time.time()
    x = x0 + d_x
    im(x,ret=False) # Only to update the mapper
    u, g2 = gm(x)
    f = -g2

    # Update hessian
    time5 = time.time()
    if update == 'powell':
        d_g = np.subtract(f0, f)
        Powell(d_x.flatten(), d_g.flatten(), h)
    elif update == 'recompute':
        get_hessian(h,gm,x)

    #Store action
    time6 = time.time()
    action[0] = gm.dforces.pot * 1/(im.temp  * im.dbeads.nbeads * units.Constants.kb)
    action[1] = im.pot / (im.temp * im.dbeads.nbeads * units.Constants.kb)
    time7 = time.time()
    # Note that for the half polymer the factor 2 cancels out (One factor in *.pot and one factor in *.nbeads)

    print ''
    print ''
    print 'Add spring %f' %(time1-time0)
    print 'Clean hessian %f' %(time2 - time1)
    print 'Nichols %f' %(time3 - time2)
    print 'Rescale %f' %(time4 - time3)
    print 'Make movement %f' %(time5 - time4)
    print 'Update hessian %f' % (time6 - time5)
    print 'Compute action %f' % (time7 - time6)
#-------------------------------------------------------------------------------------------------------------
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
            h[j + i * ii, i * ii:(i + 1) * ii] = g[i, :]

    u,g=gm(x0) #Keep the mapper updated

    ## Ask Michelle about transfer force here II
    #gm.dbeads.q = ddbeads.q
    #gm.dforces.transfer_forces(ddforces)


def spring_hessian(im):
    """Compute the 'spring hessian'           
       IN     im      = instanton mapper
       
       OUT    h       = hessian with only the spring terms ('spring hessian')
        """
    info(" @spring_hessian", verbosity.high)
    ii = im.dbeads.natoms * 3
    h = np.zeros([ii * im.dbeads.nbeads, ii * im.dbeads.nbeads])

    if im.dbeads.nbeads == 1:
        return h

    # Diagonal
    h_sp = im.dbeads.m3[0] * im.omega2
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

    return h

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
def clean_hessian(h,q, natoms,nbeads,m,m3,asr,mofi=False):
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
            w      = dynmatrix with the external modes projected out
            
        #Adapted from ipi/engine/motion/phonons.py apply_asr    """
    info(" @clean_hessian", verbosity.high)
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

           I,U=(np.linalg.eig(moi))
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
    #print dd[0:9]
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

    if nzero.size > 0:
        if np.linalg.norm(nzero) > 0.01:
            info(" Warning @Clean hessian: We have deleted %d 'zero' frequencies " % (nzero.size), verbosity.high)
            info(" but the norm is greater than 0.01 cm^-1.  This should not happen." % (nzero.size), verbosity.high)

        d = np.delete(d,range(nneg.size,nneg.size+nzero.size))
        w = np.delete(w, range(nneg.size,nneg.size+nzero.size),axis=1)

    if mofi and asr=='poly':
        return d,w,np.prod(I)
    else:
        return d,w

#----------------------------------------------------------------------------------------------------------

def get_doble_hessian(h0,im):
    """Takes a hessian of the half polymer (only the physical part) and construct the 
       hessian for the full polymer (physical + spring).
       
       IN  h0      = physical hessian of the half polymer
            im     = instanton mapper
       
       OUT h      = full ring polymer hessian (physical + spring terms)"""
    info(" @get_doble_hessian" , verbosity.high)

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


    # Spring part
    h_sp = im.dbeads.m3[0] * im.omega2

    # Diagonal
    diag = np.diag(2.0 * h_sp)

    for i in range(0, 2*im.dbeads.nbeads):
        h[i * ii:(i + 1) * ii, i * ii:(i + 1) * ii] += diag

    # quasi-band
    ndiag = np.diag(-h_sp)
    for i in range(0, 2*im.dbeads.nbeads - 1):
        h[i * ii:(i + 1) * ii, (i + 1) * ii:(i + 2) * ii] += ndiag
        h[(i + 1) * ii:(i + 2) * ii, i * ii:(i + 1) * ii] += ndiag

    # Corner
    h[0:ii, (2*nbeads - 1) * ii:(2*nbeads) * ii] += ndiag
    h[(2*nbeads - 1) * ii:(2*nbeads) * ii, 0:ii] += ndiag

    return h

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

    if im.mode == 'half':
        factor=2.0
    elif im.mode == 'full':
        factor=1.0

    betaP = 1.0 / (im.temp  * factor*im.dbeads.nbeads * units.Constants.kb)

    print >> outfile, "# Instanton calculation:"
    print >> outfile, ('  ')
    print >> outfile, ('We have %i beads in the mode  "%s" and the temperature  is %f K ' %(im.dbeads.nbeads,im.mode,units.unit_to_user('temperature',"kelvin",im.temp)))
    print >> outfile, ('Beta in a.u.', 1/im.temp)

    #action1 = gm.dforces.pot * 1/(im.temp  * im.dbeads.nbeads * units.Constants.kb)
    #action2 = im.pot / (im.temp * im.dbeads.nbeads * units.Constants.kb)
    # Note that for the half polymer the factor 2 cancels out (One factor in *.pot and one factor in *.nbeads)
    action1 = action[0]
    action2 = action[1]
    action = action1 + action2
    if im.mode == 'half':
        BN =  2*np.sum(gm.dbeads.m3[1:,:]*(gm.dbeads.q[1:,:] - gm.dbeads.q[:-1,:])**2)
    if im.mode == 'full':
        BN = np.sum(gm.dbeads.m3*(np.roll(gm.dbeads.q,1,axis=0) - gm.dbeads.q) ** 2)

    print >> outfile, ('S1/hbar   '  + str(action1 ) )
    print >> outfile, ('S2/hbar   '  + str(action2 ) )
    print >> outfile, ('S/hbar    '  + str(action  ) )
    print >> outfile, ('BN   ' + (str(BN)))
    print >> outfile, ('BN*N   ' + (str(BN * factor * im.dbeads.nbeads)))
    print >> outfile, ('  ')

    info(" ", verbosity.low)
    info(" We are almost done.", verbosity.low)
    if rates != 'true':
        info(" Compute rates is false, we are not going to compute them.", verbosity.low)
    else:
        info(" Computing rates. For this we need to compute the hessian.", verbosity.low)

        if im.mode =='full':
            get_hessian(h, gm, im.dbeads.q)
            if gm.dbeads.nbeads != 1:
                hf = np.add(im.h,h)
            else:
                hf=h.copy()
            d,w,detI = clean_hessian(hf,im.dbeads.q,im.dbeads.natoms,im.dbeads.nbeads,im.dbeads.m,im.dbeads.m3,asr,mofi=True)
            m = im.dbeads.m
        elif im.mode =='half':
            get_hessian(h, gm, im.dbeads.q)
            hbig = get_doble_hessian(h,im)
            q,nbeads,m,m3 = get_doble(im.dbeads.q, im.dbeads.nbeads, im.dbeads.m, im.dbeads.m3)
            d, w,detI = clean_hessian(hbig, q, im.dbeads.natoms, nbeads, m, m3, asr,mofi=True)
        print >> outfile, "Final  lowest ten frequencies (cm^-1)"
        print >> outfile, str( np.sign(d[0:10]) * np.absolute(d[0:10]) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17) ) # convert to cm^-1
        print >> outfile, ('  ')

        #Compute Qrot
        #if asr=='poly':
        #    qrot = ( 8*np.pi*detI / ( (units.Constants.hbar)**6 * (betaP)**3 ))**0.5
        #else:
        #    qrot = 1.0
        #print >> outfile, ('DetI ' + str(detI))
        #print >> outfile, ('Qrot ' + str(qrot))

        #Compute Qtras
        #qtras= ( ( np.sum(m)*im.dbeads.nbeads*factor ) / ( 2*np.pi*betaP*units.Constants.hbar**2 ) )**1.5
        #print >> outfile, ('Qtras ' + str(qtras))

        #Compute loqQvib
        if im.dbeads.nbeads >1: #instanton calculation
           logQvib = np.sum( np.log( betaP*units.Constants.hbar*np.sqrt(np.absolute(np.delete(d,1))) ))
           print 'Deleted frequency for computing Qvib  %f cm^-1' % (np.sign(d[1]) * np.absolute(d[1]) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17))
        else: # TS calculation,exact expression
           logQvib = np.sum( np.log(2 * np.sinh((betaP *units.Constants.hbar * np.sqrt(np.delete(d, 0)) / 2.0))))
           #logQvib = np.sum( np.log( betaP * units.Constants.hbar * np.sqrt(np.absolute(np.delete(d, 0)))))
           print 'Deleted frequency for computing Qvib  %f cm^-1' % (np.sign(d[0]) * np.absolute(d[0]) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17))

        print >> outfile, ('logQvib ' + str(logQvib))

    outfile.close()
