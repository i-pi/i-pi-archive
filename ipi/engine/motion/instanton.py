"""
Contains classes for different geometry optimization algorithms.

TODO

Algorithms implemented by Michele Ceriotti and Benjamin Helfrecht, 2015
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
from ipi.utils.instanton import Instanton, expand_hessian,add_spring_hessian,get_imvector,print_instanton,get_hessian 
from ipi.utils import units

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
                 biggest_step=1.0,
                 old_pos=np.zeros(0, float),
                 old_pot= np.zeros(0, float),
                 old_force=np.zeros(0, float),
                 hessian=np.eye(0, 0, 0, float),
                 tolerances={"energy": 1e-8, "force": 1e-8, "position": 1e-8},
                 delta=np.zeros(0, float),
                 min_temp=np.zeros(0, float),
                 hessian_init=None,
                 hessian_update=None,
                 hessian_asr=None,
                 final_rates=None):

        """Initialises InstantonMotion.
        """

        super(InstantonMotion, self).__init__(fixcom=fixcom, fixatoms=fixatoms)

        # Optimization Options
        self.mode         = mode
        self.big_step     = biggest_step
        self.tolerances   = tolerances

        self.old_x        = old_pos
        self.old_u        = old_pot
        self.old_f        = old_force

        #
        self.optimizer      = InstantonOptimizer()
        self.hessian        = hessian
        self.hessian_init   = hessian_init
        self.hessian_update = hessian_update
        self.hessian_asr    = hessian_asr
        self.final_rates    = final_rates
        self.delta          = delta
        self.temp           = min_temp


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

        if self.mode=='full':
            self.dbeads.q = x.copy()
            e=self.dbeads.vpath*self.omega2
            g=-self.dbeads.fpath*self.omega2
            self.save(e,g)

        elif self.mode=='half':
            self.dbeads.q = x.copy()
            e = np.zeros(1)
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

        if geop.temp == 0:
            if self.beads.nbeads != 1:
                raise ValueError("Temperature must be specified for an Instanton calculation ")

        self.tolerances     = geop.tolerances
        self.big_step       = geop.big_step
        self.delta          = geop.delta
        self.temp           = geop.temp

        self.hessian_update = geop.hessian_update
        self.hessian_asr    = geop.hessian_asr
        self.hessian_init   = geop.hessian_init
        self.final_rates    = geop.final_rates

        self.gm.bind(self)
        self.im.bind(self)

        if geop.hessian.size != (self.beads.q.size * self.beads.q.size):
            if geop.hessian.size == 0: #Hessian not provided
                if geop.hessian_init =='true' and self.beads.nbeads == 1:
                     info(" Initial classical hessian is not provided. We are going to compute it.", verbosity.low)
                     geop.hessian = np.zeros(self.beads.q.size, self.beads.q.size)
                else:
                    raise ValueError("nbeads >1 and/or hessian_init != 'True'  --->  An initial hessian must be provided")
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
        info("\nMD STEP %d" % step, verbosity.debug)
        info("\nMD STEP %d" % step, verbosity.low) #ALBERTO


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
                        self.beads.q[i, :] += self.delta * np.cos((i+1) *2.0*np.pi / float(self.beads.nbeads)) * imvector[:]
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
                    expand_hessian(self.hessian_init,self.hessian,self.beads.natoms)

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
   # self.exitstep(self.forces.pot, self.old_u, d_x_max)

        fx= self.forces.pot
        u0= self.old_u
        x=d_x_max
        u,g=self.gm(self.gm.dbeads.q)

    #def exitstep(self, fx, u0, x):

        #Modification of exitstep !ALBERTO spring energy
        """ Exits the simulation step. Computes time, checks for convergence. """
        info(" @GEOP: Updating bead positions", verbosity.debug)


        self.qtime += time.time()
        print "Energy", np.absolute((fx - u0) / self.beads.natoms), self.tolerances["energy"]
        print "Abs_for", np.amax(np.absolute(self.forces.f+self.im.f)), self.tolerances["force"]
        print "pos", x, self.tolerances["position"]
        if (np.absolute((fx - u0) / self.beads.natoms) <= self.tolerances["energy"]) \
                and ((np.amax(np.absolute(self.forces.f + self.im.f)) <= self.tolerances["force"]) or
                         (np.linalg.norm(self.forces.f.flatten() - self.old_f.flatten()) <= 1e-08)) \
                and (x <= self.tolerances["position"]):
            print_instanton(self.hessian, self.gm,self.im,self.hessian_asr, self.final_rates)
            softexit.trigger("Geometry optimization converged. Exiting simulation")