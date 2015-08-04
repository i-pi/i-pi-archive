"""Contains the classes that deal with the different dynamics required in
different types of ensembles.

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


Holds the algorithms required for normal mode propagators, and the objects to
do the constant temperature and pressure algorithms. Also calculates the
appropriate conserved energy quantity for the ensemble of choice.

"""

__all__=['NEBMover']

import numpy as np
import time


from ipi.engine.mover import Mover
from ipi.utils.depend import *
from ipi.utils import units
from ipi.utils.softexit import softexit
from ipi.utils.io import read_file
from ipi.utils.io.inputs.io_xml import xml_parse_file
from ipi.utils.units import Constants, unit_to_internal
from ipi.utils.mintools import min_brent, min_approx, BFGS, L_BFGS, L_BFGS_nls, min_brent_neb
from ipi.utils.messages import verbosity, warning, info

class NEBLineMover(object):
   """Creation of the one-dimensional function that will be minimized"""
   
   def __init__(self):
      self.x0 = self.d = None
      self.kappa = None
      self.first = True

   def bind(self, ens):
      self.dbeads = ens.beads.copy()
      self.dcell = ens.cell.copy()
      self.dforces = ens.forces.copy(self.dbeads, self.dcell)      
      
   def set_dir(self, x0, mdir):      
      self.x0 = x0.copy()
      self.d = mdir.copy()/np.sqrt(np.dot(mdir.flatten(),mdir.flatten()))
      if self.x0.shape != self.d.shape: raise ValueError("Incompatible shape of initial value and displacement direction")      
   
   def __call__(self, x):
      if self.first is True:
          self.dbeads.q = x
      else:
          self.dbeads.q = self.x0 + self.d * x 
      bq = depstrip(self.dbeads.q).copy()
      bf = depstrip(self.dforces.f).copy()
      be = depstrip(self.dforces.pot).copy()

      nimg = self.dbeads.nbeads
      nat = self.dbeads.natoms
      kappa = np.zeros(nimg)
      
      # get tangents, end images are distinct, fixed, pre-relaxed configurations
      btau = np.zeros((nimg, 3 * nat), float)
      for ii in range(1, nimg - 1):
          d1 = bq[ii] - bq[ii - 1] #tau plus
          d2 = bq[ii + 1] - bq[ii] #tau minus
          
      # "Old" implementation of NEB
      btau[ii] = d1 / np.linalg.norm(d1) + d2 / np.linalg.norm(d2)
      btau[ii] *= 1.0 / np.linalg.norm(btau)
          
#          # Energy of images: (ii+1) < (ii) < (ii-1)
#          if (be[ii + 1] < be[ii]) and (be[ii] < be[ii - 1]):
#              btau[ii] = d2
#
#          # Energy of images (ii-1) < (ii) < (ii+1)
#          elif (be[ii - 1] < be[ii]) and (be[ii] < be[ii + 1]):
#              btau[ii] = d1
#          
#          # Energy of image (ii) is a minimum or maximum
#          else:
#              maxpot = max(be[ii + 1] - be[ii], be[ii - 1], be[ii])
#              minpot = min(be[ii + 1] - be[ii], be[ii - 1], be[ii])
#              
#              if be[ii + 1] < be[ii - 1]:
#                  btau[ii] = d1 * minpot + d2 * maxpot
#
#              elif be[ii - 1] < be[ii + 1]:
#                  btau[ii] = d1 * maxpot + d2 * minpot
#
#              else:
#                  print "Error in NEB tangents: Energy of images are equal"
#          
#          btau[ii] *= 1.0 / np.linalg.norm(btau) 


      #if mode == "variablesprings": #TODO: input option for variable spring mode
      
          
#      if mode == "ci":
#
#      # Climbing NEB term. Choose highest energy bead after 5 (arbitrary) iterations
#          if step >= 5:
#              imax = np.argmax(be)
#              bf[imax] = bf[imax] - 2 * np.dot(bf[imax], btau[imax]) * btau[imax] 
#          
#              # Determine variable spring constants
#              #kappa = np.zeros(nimg)
#              #ei = np.zeros(nimg)
#              #emax = np.amax(be)
#              #eref = max(be[0], be[nimg])
#              #kappamax = self.kappa_max
#              #kappamin = self.kappa_min #TODO: input options for max and min spring constant
#              #deltakappa = kappamax - kappamin
#              #for ii in range(1, nimg - 1):
#              #    ei[ii] = max(be[ii], be[ii - 1])
#              #    if ei[j] > eref:
#              #        kappa[ii] = kappamax - deltakappa * ((emax - ei[ii]) / (emax - eref)) 
#              #    else:
#              #        kappa[ii] = kappamin
#              
#      else:    
#          kappa.fill(self.kappa)
#      
#       
#          # get perpendicular forces 
#          for ii in range(1, nimg - 1):
#              bf[ii] = bf[ii] - np.dot(bf[ii], btau[ii]) * btau[ii]
#          
#          # adds the spring forces           
#          for ii in range(1, nimg - 1):
#              bf[ii] += kappa[ii] * btau[ii] * np.dot(btau[ii], (bq[ii + 1] + bq[ii - 1] - 2 * bq[ii])) 

      kappa.fill(self.kappa)
       
      # get perpendicular forces 
      for ii in range(1, nimg - 1):
          bf[ii] = bf[ii] - np.dot(bf[ii], btau[ii]) * btau[ii]
          
      # adds the spring forces           
      for ii in range(1, nimg - 1):
          bf[ii] += kappa[ii] * btau[ii] * np.dot(btau[ii], (bq[ii + 1] + bq[ii - 1] - 2 * bq[ii])) 

      if self.first is True:
          self.d = bf
          self.first = False

      force = bf
      g = -np.dot(bf.flatten(), self.d.flatten())
      g = abs(g) 
      linefunc = (force, g)
      return linefunc       


class NEBBFGSMover(object):
   """ Creation of the multi-dimensional function that will be minimized"""

   def __init__(self):
      self.x0 = self.d = self.xold = None
      self.kappa = None

   def bind(self, ens):
      self.dbeads = ens.beads.copy()
      self.dcell = ens.cell.copy()
      self.dforces = ens.forces.copy(self.dbeads, self.dcell)

   def __call__(self, x):
      self.dbeads.q = x
      bq = depstrip(self.dbeads.q).copy()
      bf = depstrip(self.dforces.f).copy()
      be = depstrip(self.dforces.pot).copy()

      nimg = self.dbeads.nbeads
      nat = self.dbeads.natoms
      kappa = np.zeros(nimg)
      
      # get tangents, end images are distinct, fixed, pre-relaxed configurations
      btau = np.zeros((nimg, 3 * nat), float)
      for ii in range(1, nimg - 1):
          d1 = bq[ii] - bq[ii - 1] #tau plus
          d2 = bq[ii + 1] - bq[ii] #tau minus
          
      # "Old" implementation of NEB
      btau[ii] = d1 / np.linalg.norm(d1) + d2 / np.linalg.norm(d2)
      btau[ii] *= 1.0 / np.linalg.norm(btau)
          
#          # Energy of images: (ii+1) < (ii) < (ii-1)
#          if (be[ii + 1] < be[ii]) and (be[ii] < be[ii - 1]):
#              btau[ii] = d2
#
#          # Energy of images (ii-1) < (ii) < (ii+1)
#          elif (be[ii - 1] < be[ii]) and (be[ii] < be[ii + 1]):
#              btau[ii] = d1
#          
#          # Energy of image (ii) is a minimum or maximum
#          else:
#              maxpot = max(be[ii + 1] - be[ii], be[ii - 1], be[ii])
#              minpot = min(be[ii + 1] - be[ii], be[ii - 1], be[ii])
#              
#              if be[ii + 1] < be[ii - 1]:
#                  btau[ii] = d1 * minpot + d2 * maxpot
#
#              elif be[ii - 1] < be[ii + 1]:
#                  btau[ii] = d1 * maxpot + d2 * minpot
#
#              else:
#                  print "Error in NEB tangents: Energy of images are equal"
#          
#          btau[ii] *= 1.0 / np.linalg.norm(btau) 


      #if mode == "variablesprings": #TODO: input option for variable spring mode
      
          
#      if mode == "ci":
#
#      # Climbing NEB term. Choose highest energy bead after 5 (arbitrary) iterations
#          if step >= 5:
#              imax = np.argmax(be)
#              bf[imax] = bf[imax] - 2 * np.dot(bf[imax], btau[imax]) * btau[imax] 
#          
#              # Determine variable spring constants
#              #kappa = np.zeros(nimg)
#              #ei = np.zeros(nimg)
#              #emax = np.amax(be)
#              #eref = max(be[0], be[nimg])
#              #kappamax = self.kappa_max
#              #kappamin = self.kappa_min #TODO: input options for max and min spring constant
#              #deltakappa = kappamax - kappamin
#              #for ii in range(1, nimg - 1):
#              #    ei[ii] = max(be[ii], be[ii - 1])
#              #    if ei[j] > eref:
#              #        kappa[ii] = kappamax - deltakappa * ((emax - ei[ii]) / (emax - eref)) 
#              #    else:
#              #        kappa[ii] = kappamin
#              
#      else:    
#          kappa.fill(self.neb_kappa)
#      
      kappa.fill(self.kappa)

      # get perpendicular forces 
      for ii in range(1, nimg - 1):
          bf[ii] = bf[ii] - np.dot(bf[ii], btau[ii]) * btau[ii]
          
      # adds the spring forces           
      for ii in range(1, nimg - 1):
          bf[ii] += kappa[ii] * btau[ii] * np.dot(btau[ii], (bq[ii + 1] + bq[ii - 1] - 2 * bq[ii])) 

      e = self.dforces.pot #0.0
      g = -bf
      return e, g        

class NEBMover(Mover):
   """Nudged elastic band routine.

   Attributes:

   """

   def __init__(self, fixcom=False, fixatoms=None,
             mode="sd", 
             grad_tolerance=1.0e-6, maximum_step=100.0,
             cg_old_force=np.zeros(0, float),
             cg_old_direction=np.zeros(0, float),
             invhessian=np.eye(0), 
             ls_options={ "tolerance": 1e-5,  "iter": 100.0 , "step": 1e-3, "adaptive":1.0 } ,
             tolerances={"energy" : 1e-5, "force": 1e-5, "position": 1e-5},
             corrections = 5,
             qlist = np.zeros(0, float),
             glist = np.zeros(0, float),
             interior = False,
             endpoints = True,
             spring = {"varsprings": False, "kappa": 1.0, "kappamax": 1.5, "kappamin": 0.5},
             climb = False
             ):   
                 
      """Initialises NEBMover.

      Args:
         fixcom: An optional boolean which decides whether the centre of mass
            motion will be constrained or not. Defaults to False.         
      """

      super(NEBMover,self).__init__(fixcom=fixcom, fixatoms=fixatoms)
      
      # optimization options
      self.ls_options = ls_options
      self.tolerances = tolerances
      self.mode = mode
      self.grad_tol = grad_tolerance
      self.max_step = maximum_step
      self.cg_old_f = cg_old_force
      self.cg_old_d = cg_old_direction
      self.invhessian = invhessian
      self.corrections = corrections
      self.qlist = qlist
      self.glist = glist
      self.interior = interior
      self.endpoints = endpoints
      self.spring = spring # TODO: should contain all spring constants
      self.climb = climb

      self.neblm = NEBLineMover()
      self.nebbfgsm = NEBBFGSMover()
   
   def bind(self, beads, nm, cell, bforce, bbias, prng):
      
      super(NEBMover,self).bind(beads, nm, cell, bforce, bbias, prng)
      if self.cg_old_f.shape != beads.q.shape :
         if self.cg_old_f.shape == (0,): 
            self.cg_old_f = np.zeros(beads.q.shape, float)
         else: 
            raise ValueError("Conjugate gradient force size does not match system size")
      if self.cg_old_d.shape != beads.q.shape :
         if self.cg_old_d.shape == (0,): 
            self.cg_old_d = np.zeros(beads.q.shape, float)
         else: 
            raise ValueError("Conjugate gradient direction size does not match system size")
      if self.invhessian.size != (beads.q.size * beads.q.size):
          if self.invhessian.size == 0:
              self.invhessian = np.eye(beads.q.size, beads.q.size, 0, float)
          else:
              raise ValueError("Inverse Hessian size does not match system size")

      self.neblm.bind(self)
      self.nebbfgsm.bind(self)
            
   def step(self, step=None):
      """Does one simulation time step."""

      info("\nMD STEP %d" % step, verbosity.debug)
      self.nebbfgsm.kappa = self.spring["kappa"]
      self.neblm.kappa = self.spring["kappa"]

      self.ptime = self.ttime = 0
      self.qtime = -time.time()

      if self.mode == "lbfgs":

          # L-BFGS Minimization
          # Initialize approximate Hessian inverse to the identity, and direction
          # to the steepest descent direction
          # Initialize lists of previous positions and gradients
          
          if step == 0: # or np.sqrt(np.dot(self.bfgsm.d, self.bfgsm.d)) == 0.0: <-- this part for restarting at claimed minimum
              info(" @GEOP: Initializing L-BFGS", verbosity.debug)
              fx, nebgrad = self.nebbfgsm(self.beads.q) #depstrip(self.forces.f) / np.sqrt(np.dot(self.forces.f.flatten(), self.forces.f.flatten()))
              self.nebbfgsm.d = -nebgrad 
              self.nebbfgsm.xold = self.beads.q.copy()
              self.qlist = np.zeros((self.corrections, len(self.beads.q.flatten())))
              self.glist = np.zeros((self.corrections, len(self.beads.q.flatten())))

          # Current energy, forces, and function definitions
          # for use in finding converged minimization

          else:
              fx, nebgrad = self.nebbfgsm(self.beads.q) #depstrip(self.forces.f) / np.sqrt(np.dot(self.forces.f.flatten(), self.forces.f.flatten()))
          u0, du0 = (fx, nebgrad)
          self.cg_old_f[:] = -nebgrad
          
          # If we have not yet stored all requested corrections, do normal BFGS and store
          # the corrections
          self.beads.q, fx, self.nebbfgsm.d, self.qlist, self.glist = L_BFGS_nls(self.beads.q, 
                self.nebbfgsm.d, self.nebbfgsm, self.qlist, self.glist, 
                fdf0=(u0, du0), max_step=self.max_step, tol=self.ls_options["tolerance"], 
                grad_tol=self.ls_options["gradtolerance"], itmax=self.ls_options["iter"], 
                init_step=self.ls_options["step"],
                m=self.corrections, k=step)

          info(" @GEOP: Updated position list", verbosity.debug)
          info(" @GEOP: Updated gradient list", verbosity.debug)

          x = np.amax(np.absolute(np.subtract(self.beads.q, self.nebbfgsm.xold)))
          self.nebbfgsm.xold[:] = self.beads.q

          info(" @GEOP: Updated bead positions", verbosity.debug)

      # Routine for steepest descent and conjugate gradient
      else:
          if (self.mode == "sd" or step == 0): 
          
              # Steepest descent minimization
              # gradf1 = force at current atom position
              # dq1 = direction of steepest descent
              # dq1_unit = unit vector of dq1
              nebgrad = self.neblm(self.beads.q)[0]
              gradf1 = dq1 = -nebgrad

              # move direction for steepest descent and 1st conjugate gradient step
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
              nebgrad = self.neblm(self.beads.q)[0]
              gradf1 = -nebgrad 
              beta = np.dot((gradf1.flatten() - gradf0.flatten()), gradf1.flatten()) / (np.dot(gradf0.flatten(), gradf0.flatten()))
              dq1 = gradf1 + max(0.0, beta) * dq0
              dq1_unit = dq1 / np.sqrt(np.dot(dq1.flatten(), dq1.flatten()))
              info(" @GEOP: Determined CG direction", verbosity.debug)

          # store force and direction for next CG step1
          self.cg_old_d[:] = dq1
          self.cg_old_f[:] = gradf1
   
          if (len(self.fixatoms)>0):
              for dqb in dq1_unit:
                  dqb[self.fixatoms*3] = 0.0
                  dqb[self.fixatoms*3+1] = 0.0
                  dqb[self.fixatoms*3+2] = 0.0
      
          self.neblm.set_dir(depstrip(self.beads.q), dq1_unit)

          # reuse initial value since we have energy and forces already
          u0 = np.dot(-nebgrad.flatten(), dq1_unit.flatten())
          u0 = np.sqrt(np.dot(u0, u0))

          (x, nebgrad_min) = min_brent_neb(self.neblm, fdf0=u0, x0=0.0, 
                  tol=self.ls_options["tolerance"], 
                  itmax=self.ls_options["iter"], init_step=self.ls_options["step"]) 

          # automatically adapt the search step for the next iteration. 
          # relaxes better with very small step --> multiply by factor of 0.1 or 0.01
          self.ls_options["step"] = 0.1 * x * self.ls_options["adaptive"] + (1 - self.ls_options["adaptive"]) * self.ls_options["step"] 
      
          self.beads.q += dq1_unit * x
          info(" @GEOP: Updated bead positions", verbosity.debug)

      self.qtime += time.time()
      
#      # Determine conditions for converged relaxation
#      if ((fx - u0) / self.beads.natoms <= self.tolerances["energy"])\
#          and ((np.amax(np.absolute(self.forces.f)) <= self.tolerances["force"])\
#              or (np.sqrt(np.dot(self.forces.f.flatten() - self.cg_old_f.flatten(),\
#                  self.forces.f.flatten() - self.cg_old_f.flatten())) == 0.0))\
#          and (x <= self.tolerances["position"]):
#          softexit.trigger("Geometry optimization converged. Exiting simulation")
