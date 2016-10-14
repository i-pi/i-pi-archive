"""
TODO

Algorithms implemented by Michele Ceriotti and Benjamin Helfrecht, 2015
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np
import time

from ipi.engine.motion import Motion
from ipi.utils.depend import depstrip
from ipi.utils.softexit import softexit
from ipi.utils.mintools import min_brent, BFGS, L_BFGS
from ipi.utils.messages import verbosity, info


__all__ = ['GeopMover']


class LineMover(object):
    """Creation of the one-dimensional function that will be minimized.
    Used in steepest descent and conjugate gradient minimizers.

    Attributes:
        x0: initial position
        d: move direction
    """

    def __init__(self):
        self.x0 = self.d = None

    def bind(self, ens):
        self.dbeads = ens.beads.copy()
        self.dcell = ens.cell.copy()
        self.dforces = ens.forces.copy(self.dbeads, self.dcell)

    def set_dir(self, x0, mdir):
        self.x0 = x0.copy()
        self.d = mdir.copy() / np.sqrt(np.dot(mdir.flatten(), mdir.flatten()))
        if self.x0.shape != self.d.shape:
            raise ValueError("Incompatible shape of initial value and displacement direction")

    def __call__(self, x):

        self.dbeads.q = self.x0 + self.d * x
        e = self.dforces.pot   # Energy
        g = - np.dot(depstrip(self.dforces.f).flatten(), self.d.flatten())   # Gradient
        return e, g


class BFGSMover(object):
    """Creation of the multi-dimensional function that will be minimized.
    Used in the BFGS and L-BFGS minimizers.

    Attributes:
        x0: initial position
        d: move direction
        xold: previous position
    """

    def __init__(self):
        self.x0 = None
        self.d = None
        self.xold = None

    def bind(self, ens):
        self.dbeads = ens.beads.copy()
        self.dcell = ens.cell.copy()
        self.dforces = ens.forces.copy(self.dbeads, self.dcell)

    def __call__(self, x):
        self.dbeads.q = x
        e = self.dforces.pot   # Energy
        g = - self.dforces.f   # Gradient
        return e, g


class GeopMover(Motion):
    """Geometry optimization routine. Controls direction choice for
    steepest descent and conjugate gradient. Checks for satisfaction of
    tolerances to exit minimization.

    Attributes:
        mode: minimization algorithm to use
        maximum_step: max allowed step size for BFGS/L-BFGS
        cg_old_force: force on previous step
        cg_old_direction: move direction on previous step
        invhessian: stored inverse Hessian matrix for BFGS
        ls_options:
            tolerance: energy tolerance for exiting minimization algorithm
            iter: maximum number of allowed iterations for minimization algorithm for each MD step
            step: initial step size for steepest descent and conjugate gradient
            adaptive: T/F adaptive step size for steepest descent and conjugate
                gradient
        tolerances:
            energy: change in energy tolerance for ending minimization
            force: force/change in force tolerance foe ending minimization
            position: change in position tolerance for ending minimization
        corrections: number of corrections to be stored for L-BFGS
        qlist: list of previous positions (x_n+1 - x_n) for L-BFGS. Number of entries = corrections
        glist: list of previous gradients (g_n+1 - g_n) for L-BFGS. Number of entries = corrections
    """

    def __init__(self, fixcom=False, fixatoms=None,
                 mode="sd",
                 maximum_step=100.0,
                 cg_old_force=np.zeros(0, float),
                 cg_old_direction=np.zeros(0, float),
                 invhessian=np.eye(0, 0, 0, float),
                 ls_options={"tolerance": 1e-6, "iter": 100.0, "step": 1e-3, "adaptive": 1.0},
                 tolerances={"energy": 1e-8, "force": 1e-8, "position": 1e-8},
                 corrections=5,
                 qlist=np.zeros(0, float),
                 glist=np.zeros(0, float)):
        """Initialises GeopMover.

        Args:
           fixcom: An optional boolean which decides whether the centre of mass
              motion will be constrained or not. Defaults to False.
        """

        super(GeopMover, self).__init__(fixcom=fixcom, fixatoms=fixatoms)

        # optimization options
        self.ls_options = ls_options
        self.tolerances = tolerances
        self.mode = mode
        self.max_step = maximum_step
        self.cg_old_f = cg_old_force
        self.cg_old_d = cg_old_direction
        self.invhessian = invhessian
        self.corrections = corrections
        self.qlist = qlist
        self.glist = glist

        self.lm = LineMover()
        self.bfgsm = BFGSMover()

    def bind(self, ens, beads, nm, cell, bforce, prng):

        super(GeopMover, self).bind(ens, beads, nm, cell, bforce, prng)
        if self.cg_old_f.shape != beads.q.size:
            if self.cg_old_f.size == 0:
                self.cg_old_f = np.zeros(beads.q.size, float)
            else:
                raise ValueError("Conjugate gradient force size does not match system size")
        if self.cg_old_d.size != beads.q.size:
            if self.cg_old_d.size == 0:
                self.cg_old_d = np.zeros(beads.q.size, float)
            else:
                raise ValueError("Conjugate gradient direction size does not match system size")
        if self.invhessian.size != (beads.q.size * beads.q.size):
            if self.invhessian.size == 0:
                self.invhessian = np.eye(beads.q.size, beads.q.size, 0, float)
            else:
                raise ValueError("Inverse Hessian size does not match system size")

        self.lm.bind(self)
        self.bfgsm.bind(self)

    def step(self, step=None):
        """Does one simulation time step."""

        self.ptime = 0.0
        self.ttime = 0.0
        self.qtime = -time.time()

        info("\nMD STEP %d" % step, verbosity.debug)

        if self.mode == "bfgs":

            # BFGS Minimization
            # Initialize approximate Hessian inverse to the identity and direction
            # to the steepest descent direction
            if step == 0:   # or np.sqrt(np.dot(self.bfgsm.d, self.bfgsm.d)) == 0.0: <-- this part for restarting at claimed minimum (optional)
                info(" @GEOP: Initializing BFGS", verbosity.debug)
                self.bfgsm.d = depstrip(self.forces.f) / np.sqrt(np.dot(self.forces.f.flatten(), self.forces.f.flatten()))
                self.bfgsm.xold = self.beads.q.copy()

            # Current energy and forces
            u0 = self.forces.pot.copy()
            du0 = - self.forces.f

            # Store previous forces
            self.cg_old_f[:] = self.forces.f

            # Do one iteration of BFGS, return new point, function value,
            # move direction, and current Hessian to use for next iteration
            self.beads.q, fx, self.bfgsm.d, self.invhessian = BFGS(self.beads.q,
                self.bfgsm.d, self.bfgsm, fdf0=(u0, du0), invhessian=self.invhessian,
                max_step=self.max_step, tol=self.ls_options["tolerance"],
                itmax=self.ls_options["iter"])

            # x = current position - previous position; use for exit tolerance
            x = np.amax(np.absolute(np.subtract(self.beads.q, self.bfgsm.xold)))

            # Store old position
            self.bfgsm.xold[:] = self.beads.q

            info(" @GEOP: Updating bead positions", verbosity.debug)

        elif self.mode == "lbfgs":

            # L-BFGS Minimization
            # Initialize approximate Hessian inverse to the identity and direction
            # to the steepest descent direction
            # Initialize lists of previous positions and gradient
            if step == 0:   # or np.sqrt(np.dot(self.bfgsm.d, self.bfgsm.d)) == 0.0: <-- this part for restarting at claimed minimum (optional)
                info(" @GEOP: Initializing L-BFGS", verbosity.debug)
                self.bfgsm.d = depstrip(self.forces.f) / np.sqrt(np.dot(self.forces.f.flatten(), self.forces.f.flatten()))
                self.bfgsm.xold = self.beads.q.copy()
                self.qlist = np.zeros((self.corrections, len(self.beads.q.flatten())))
                self.glist = np.zeros((self.corrections, len(self.beads.q.flatten())))

            # Current energy and force
            u0, du0 = (self.forces.pot.copy(), - self.forces.f)

            # Store previous forces
            self.cg_old_f[:] = self.forces.f

            # Do one iteration of L-BFGS, return new point, function value,
            # move direction, and current Hessian to use for next iteration
            self.beads.q, fx, self.bfgsm.d, self.qlist, self.glist = L_BFGS(self.beads.q,
                self.bfgsm.d, self.bfgsm, self.qlist, self.glist,
                fdf0=(u0, du0), max_step=self.max_step, tol=self.ls_options["tolerance"],
                itmax=self.ls_options["iter"],
                m=self.corrections, k=step)

            info(" @GEOP: Updated position list", verbosity.debug)
            info(" @GEOP: Updated gradient list", verbosity.debug)

            # x = current position - old position. Used for convergence tolerance
            x = np.amax(np.absolute(np.subtract(self.beads.q, self.bfgsm.xold)))

            # Store old position
            self.bfgsm.xold[:] = self.beads.q

            info(" @GEOP: Updated bead positions", verbosity.debug)

        # Routine for steepest descent and conjugate gradient
        else:
            if (self.mode == "sd" or step == 0):

                # Steepest descent minimization
                # gradf1 = force at current atom position
                # dq1 = direction of steepest descent
                # dq1_unit = unit vector of dq1
                gradf1 = dq1 = depstrip(self.forces.f)

                # Move direction for steepest descent and 1st conjugate gradient step
                dq1_unit = dq1 / np.sqrt(np.dot(gradf1.flatten(), gradf1.flatten()))
                info(" @GEOP: Determined SD direction", verbosity.debug)

            else:

                # Conjugate gradient, Polak-Ribiere
                # gradf1: force at current atom position
                # gradf0: force at previous atom position
                # dq1 = direction to move
                # dq0 = previous direction
                # dq1_unit = unit vector of dq1
                gradf0 = self.cg_old_f
                dq0 = self.cg_old_d
                gradf1 = depstrip(self.forces.f)
                beta = np.dot((gradf1.flatten() - gradf0.flatten()), gradf1.flatten()) / (np.dot(gradf0.flatten(), gradf0.flatten()))
                dq1 = gradf1 + max(0.0, beta) * dq0
                dq1_unit = dq1 / np.sqrt(np.dot(dq1.flatten(), dq1.flatten()))
                info(" @GEOP: Determined CG direction", verbosity.debug)

            # Store force and direction for next CG step
            self.cg_old_d[:] = dq1
            self.cg_old_f[:] = gradf1

            if len(self.fixatoms) > 0:
                for dqb in dq1_unit:
                    dqb[self.fixatoms*3] = 0.0
                    dqb[self.fixatoms*3+1] = 0.0
                    dqb[self.fixatoms*3+2] = 0.0

            self.lm.set_dir(depstrip(self.beads.q), dq1_unit)

            # Reuse initial value since we have energy and forces already
            u0, du0 = (self.forces.pot.copy(), np.dot(depstrip(self.forces.f.flatten()), dq1_unit.flatten()))

            # Do one SD/CG iteration; return positions and energy
            (x, fx) = min_brent(self.lm, fdf0=(u0, du0), x0=0.0,
                    tol=self.ls_options["tolerance"],
                    itmax=self.ls_options["iter"], init_step=self.ls_options["step"])

            # Automatically adapt the search step for the next iteration.
            # Relaxes better with very small step --> multiply by factor of 0.1 or 0.01
            self.ls_options["step"] = 0.1 * x * self.ls_options["adaptive"] + (1 - self.ls_options["adaptive"]) * self.ls_options["step"]

            self.beads.q += dq1_unit * x
            info(" @GEOP: Updated bead positions", verbosity.debug)

        self.qtime += time.time()

        # Determine conditions for converged relaxation
        if ((fx - u0) / self.beads.natoms <= self.tolerances["energy"])\
                and ((np.amax(np.absolute(self.forces.f)) <= self.tolerances["force"])
                    or (np.sqrt(np.dot(self.forces.f.flatten() - self.cg_old_f.flatten(),
                        self.forces.f.flatten() - self.cg_old_f.flatten())) == 0.0))\
                and (x <= self.tolerances["position"]):
            softexit.trigger("Geometry optimization converged. Exiting simulation")
