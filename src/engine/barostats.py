"""Contains the classes that deal with constant pressure dynamics.

Contains the algorithms which propagate the position and momenta steps in the
constant pressure and constant stress ensembles. Holds the properties directly
related to these ensembles, such as the internal and external pressure and 
stress and the strain energy.

Classes:
   Barostat: Base barostat class with the generic methods and attributes.
   BaroFlexi: Deals with flexible cell dynamics. Used for NST ensembles.
   BaroRigid: Deals with rigid cell dynamics. Used for NVT and NPT ensembles.
"""

__all__ = ['Barostat', 'BaroFlexi', 'BaroRigid']

import math, time
import numpy as np
from utils.depend import *
from utils.units import *
from utils.mathtools import eigensystem_ut3x3, invert_ut3x3, exp_ut3x3, det_ut3x3
from engine.thermostats import Thermostat
from inputs.thermostats import InputThermo

class Barostat(dobject): 
   """Base barostat class.

   Gives the standard methods and attributes needed in all the barostat classes.

   Attributes:
      thermostat: A thermostat object used to keep the cell momenta at a 
         specified kinetic temperature.
      beads: A beads object giving the atoms positions
      cell: A cell object giving the system box.
      forces: A forces object giving the virial and the forces acting on 
         each bead.

   Depend objects:
      sext: The external stress tensor.
      pext: The external pressure.
      dt: The time step used in the algorithms. Depends on the simulation dt.
      temp: The simulation temperature. Higher than the system temperature by
         a factor of the number of beads. Depends on the simulation temp.
      pot: The elastic strain potential for the cell. Depends on sext, the 
         reference cell volume, and the strain.
      piext: The accumulated stress compared to the reference cell. Depends
         on the reference cell volume and vector matrix, the cell volume and 
         vector matrix, the strain and sext.
      stress: The internal stress. Depends on the cell kinetic stress tensor
         and volume and the forces virial.
      press: The internal pressure. Depends on the stress.
   """

   def __init__(self, pext=0.0, sext=None, dt=None, temp=None, thermostat=None):
      """Initialises base barostat class.

      Note that the external stress and the external pressure are synchronized.
      This makes most sense going from the stress to the pressure, but if you 
      must go in the other direction the stress is assumed to be isotropic.

      Args:
         pext: Optional float giving the external pressure. Defaults to 
            Tr(sext)/3.0
         sext: Optional array givin the external stress tensor. Defaults to 
            pext*I, where I is a 3*3 identity matrix.
         dt: Optional float giving the time step for the algorithms. Defaults
            to the simulation dt.
         temp: Optional float giving the temperature for the thermostat. 
            Defaults to the simulation temp.
         thermostat: Optional thermostat object. Defaults to Thermostat().
      """
     
      dset(self,"pext",depend_value(name="pext", value=pext))
      dset(self,"sext",depend_value(name="sext", value=pext))
      
      if thermostat is None:
         thermostat = Thermostat()
      self.thermostat = thermostat   
     
      dset(self,"dt",depend_value(name='dt'))
      if dt is None:
         self.dt = 2.0*self.thermostat.dt
      else:
         self.dt = dt
      dset(self.thermostat,"dt",   
         depend_value(name="dt", func=self.get_halfdt,
            dependencies=[dget(self,"dt")],
               dependants=dget(self.thermostat,"dt")._dependants))
           
      dset(self, "temp", depend_value(name="temp", value=temp))
      deppipe(self, "temp", self.thermostat,"temp")
      if not temp is None:
         self.temp = temp
      
   def get_halfdt(self):
      """Returns half the simulation timestep."""

      return self.dt*0.5
         
   def bind(self, beads, cell, forces):
      """Binds beads, cell and forces to the barostat.

      This takes a beads object, a cell object and a forcefield object and 
      makes them members of the barostat. It also then creates the objects that
      will hold the data needed in the barostat algorithms and the dependency 
      network.

      Args:
         beads: The beads object from which the bead positions are taken.
         cell: The cell object from which the system box is taken.
         forces: The forcefield object from which the force and virial are
            taken.
      """

      self.beads=beads
      self.cell=cell
      self.forces=forces

      dset(self,"pot",
         depend_value(name='pot', func=self.get_pot, 
            dependencies=[ dget(cell,"V0"), dget(cell,"strain"), dget(self,"sext") ]))            
      dset(self,"piext",
         depend_value(name='piext', func=self.get_piext, 
            dependencies=[ dget(cell,"V0"), dget(cell,"V"), dget(cell,"h"), dget(cell,"ih0"), dget(cell,"strain"), dget(self,"sext") ]))     
      dset(self,"kstress",
         depend_value(name='kstress', func=self.get_kstress, 
            dependencies=[ dget(beads,"q"), dget(beads,"qc"), dget(self,"temp") , dget(forces,"f") ]))
      dset(self,"stress",
         depend_value(name='stress', func=self.get_stress, 
            dependencies=[ dget(self,"kstress"), dget(cell,"V"), dget(forces,"vir") ]))
      dset(self,"press",
         depend_value(name='press', func=self.get_press, 
            dependencies=[ dget(self,"stress") ]))
      
   def pstep(self):
      """Dummy momenta propagator step."""

      pass

   def qcstep(self):
      """Dummy centroid position propagator step."""

      pass   
      
   def get_pot(self):
      """Calculates the elastic strain energy of the cell."""

      return self.cell.V0*np.trace(np.dot(self.sext, self.cell.strain))

   def get_piext(self):
      """Calculates the accumulated external stress tensor.

      This tensor is calculated with respect to the reference cell, and so
      gives a measure of the stress due to deviation from the cell for zero
      external pressure.
      """

      root = np.dot(depstrip(self.cell.h), depstrip(self.cell.ih0))
      pi = np.dot(root, depstrip(self.sext))
      
      pi = np.dot(pi, np.transpose(root))
      pi *= self.cell.V0/self.cell.V
      return pi
      
   def get_kstress(self):
      """Calculates the quantum centroid virial kinetic stress tensor 
      estimator.
      """

      kst = np.zeros((3,3),float)
      q = depstrip(self.beads.q)
      qc = depstrip(self.beads.qc)
      na3 = 3*self.beads.natoms
      for b in range(self.beads.nbeads):
         for i in range(3):
            for j in range(i,3):
               kst[i,j] -= np.dot(q[b,i:na3:3] - qc[i:na3:3], 
                  depstrip(self.forces.f[b])[j:na3:3])

      for i in range(3):
         kst[i,i] += Constants.kb*self.temp*(self.beads.natoms)
      kst *= 1.0/self.beads.nbeads
      return kst

   def get_stress(self):
      """Calculates the internal stress tensor."""
      
      return (self.kstress + self.forces.vir/float(self.beads.nbeads))/self.cell.V

   def get_press(self):
      """Calculates the internal pressure."""

      return np.trace(self.stress)/3.0


class BaroFlexi(Barostat):
   """Barostat object for flexible cell simulations.

   Propagates the relevant equations of motion to give a constant stress 
   ensemble assuming an upper triangular lattice vector matrix 
   (see P. Raiteri, J. Gale and G. Bussi, J. Phys.: Condens. Matter 23, 334213 
   (2011)). Note that the volume fluctuations are assumed to 
   only affect the centroid normal mode, and so the other modes are just 
   propagated in the same way as for constant volume ensembles.
   """

   def bind(self, beads, cell, forces):
      """Binds beads, cell and forces to the barostat.

      As for the base class bind, except the cell is also bound to the 
      thermostat.

      Args:
         beads: The beads object from which the bead positions are taken.
         cell: The cell object from which the system box is taken.
         forces: The forcefield object from which the force and virial are
            taken.
      """

      super(BaroFlexi,self).bind(beads, cell, forces)
      self.thermostat.bind(cell=self.cell)

   def pstep(self):
      """Propagates the cell and centroid momenta.

      Updates the centroid momenta as for the velocity verlet momentum step, but
      also updates the cell momenta based on the mismatch of the internal and
      external stress tensors and the motion of the centroids. Note that the 
      outer product terms are done by splitting the sum up into several
      dot product calculations to make the calculation faster.
      """

      dthalf = self.dt*0.5
      dthalf2 = dthalf**2/2.0
      dthalf3 = dthalf**3/3.0     

      L = np.zeros((3,3))
      for i in range(3):
         L[i,i] = 3.0 - i
      
      self.cell.p += dthalf*(self.cell.V*(self.stress - self.piext) + 2.0*Constants.kb*self.thermostat.temp*L)       

      m = depstrip(self.beads.m)

      fc = depstrip(self.forces.fnm[0])/math.sqrt(self.beads.nbeads)
      fx = fc[0:3*self.beads.natoms:3]
      fy = fc[1:3*self.beads.natoms:3]
      fz = fc[2:3*self.beads.natoms:3]
      fxm = fx/m
      fym = fy/m
      fzm = fz/m

      pc = depstrip(self.beads.pc)
      px = pc[0:3*self.beads.natoms:3]
      py = pc[1:3*self.beads.natoms:3]
      pz = pc[2:3*self.beads.natoms:3]
      
      cp = np.zeros((3,3),float)
      cp[0,0] = dthalf2*2.0*np.dot(fxm,px) + dthalf3*np.dot(fx,fxm)
      cp[1,1] = dthalf2*2.0*np.dot(fym,py) + dthalf3*np.dot(fy,fym)
      cp[2,2] = dthalf2*2.0*np.dot(fzm,pz) + dthalf3*np.dot(fz,fzm)
      cp[0,1] = dthalf2*(np.dot(fxm,py) + np.dot(px,fym)) + dthalf3*np.dot(fx,fym)
      cp[0,2] = dthalf2*(np.dot(fxm,pz) + np.dot(px,fzm)) + dthalf3*np.dot(fx,fzm)
      cp[1,2] = dthalf2*(np.dot(fym,pz) + np.dot(py,fzm)) + dthalf3*np.dot(fy,fzm)            
      self.cell.p += cp
      self.beads.p += self.forces.f*dthalf      
      
   def qcstep(self):
      """Propagates the cell and centroid position and momenta.

      Updates the centroid and cell momenta and postions due to the motion 
      of the cell box. Note that the dot product terms in the algorithms
      for the different centroids are done in one step by first unflattening
      the position, momenta into the x, y and z components and using matrix
      multiplication.
      """

      vel_mat = depstrip(self.cell.p)/self.cell.m

      dist_mat = vel_mat*self.dt
      exp_mat = exp_ut3x3(dist_mat)
      neg_exp_mat = invert_ut3x3(exp_mat)
      sinh_mat = 0.5*(exp_mat - neg_exp_mat)
      ips_mat = np.dot( sinh_mat, invert_ut3x3(vel_mat) )

      nat = self.beads.natoms

      pc = depstrip(self.beads.pc).reshape((nat,3)) 
      qc = depstrip(self.beads.qc).reshape((nat,3))
      m = depstrip(self.beads.m)
      m3 = np.zeros((nat,3))
      for i in range(3):
         m3[:,i] = m

      qc = np.dot(qc,exp_mat.T)+np.dot(pc/m3,ips_mat.T)
      pc = np.dot(pc,neg_exp_mat.T)

      self.beads.qnm[0,:] = qc.reshape(3*nat)*math.sqrt(self.beads.nbeads)
      self.beads.pnm[0,:] = pc.reshape(3*nat)*math.sqrt(self.beads.nbeads)
                    
      self.cell.h = np.dot(exp_mat, self.cell.h)

      
class BaroRigid(Barostat):
   """Barostat object for rigid cell simulations.

   Note that the volume fluctuations are assumed to 
   only affect the centroid normal mode, and so the other modes are just 
   propagated in the same way as for constant volume ensembles. The potential
   is calculated in a simpler way than for flexible dynamics, as there is no
   contribution from change in system box shape.
   """

   def get_pot(self):
      """Calculates the elastic strain energy of the cell."""

      return self.cell.V*self.pext
      
   def bind(self, beads, cell, forces):
      """Binds beads, cell and forces to the barostat.

      As for the base class bind, except the cell is also bound to the 
      thermostat, and as the potential has been redifined it's
      dependencies must also be updated.

      Args:
         beads: The beads object from which the bead positions are taken.
         cell: The cell object from which the system box is taken.
         forces: The forcefield object from which the force and virial are
            taken.
      """

      super(BaroRigid,self).bind(beads, cell, forces)
      self.thermostat.bind(pm=(self.cell.P, self.cell.M))
      dset(self,"pot",depend_value(name='pot', func=self.get_pot, 
          dependencies=[ dget(self.cell,"V"), dget(self,"pext")  ] ) )
      
   def pstep(self):
      """Propagates the cell and centroid momenta.

      Updates the centroid momenta as for the velocity verlet momentum step, but
      also updates the cell momenta based on the mismatch of the internal and
      external pressure and the motion of the centroids.
      """
      
      dthalf = self.dt*0.5
      dthalf2 = dthalf**2/2.0
      dthalf3 = dthalf**3/3.0     
      
      self.cell.P += dthalf*3.0*(self.cell.V*(self.press - self.pext) + 2.0*Constants.kb*self.temp)

      fc = depstrip(self.forces.fnm)[0,:]/math.sqrt(self.beads.nbeads)
      m = depstrip(self.beads.centroid.m3)
      pc = depstrip(self.beads.pc)
            
      self.cell.P += dthalf2*np.dot(pc,fc/m) + dthalf3*np.dot(fc,fc/m)
      #Should dthalf2 be dthalf2*2?
   
      self.beads.p += depstrip(self.forces.f)*dthalf      
           
   def qcstep(self):
      """Propagates the cell and centroid position and momenta.

      Updates the centroid and cell momenta and postions due to the motion 
      of the cell box.
      """

      vel = self.cell.P/self.cell.M
      exp, neg_exp = (math.exp(vel*self.dt), math.exp(-vel*self.dt))
      sinh = 0.5*(exp - neg_exp)

      pc = depstrip(self.beads.pc)
      qc = depstrip(self.beads.qc)
      m = depstrip(self.beads.centroid.m3)      
      qc *= exp
      qc += (sinh/vel)*pc/m
      pc *= neg_exp

      self.beads.qnm[0,:] = qc*math.sqrt(self.beads.nbeads)
      self.beads.pnm[0,:] = pc*math.sqrt(self.beads.nbeads)

      self.cell.V *= exp**3
