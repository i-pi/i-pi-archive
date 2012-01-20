"""Contains the classes which deal with the system box. 

Used for implementing the minimum image convention, as well as holding the 
dynamical variables and the methods required for variable cell NST and NPT 
dynamics.

Classes: 
   Cell: Base cell class with the generic methouds and attributes.
   CellFlexi: Cell class with methods for flexible cell dynamics.
   CellRigid: Cell class with methods for isotropic cell dynamics.
   RestartCell: Deals with creating the Cell object from a file, and writing
      the checkpoints.
"""

__all__ = ['Cell', 'CellFlexi', 'CellRigid', 'RestartCell']

import numpy as np
import math
from utils.depend import *
from utils.restart import *
import utils.io.io_pdb
from utils.mathtools import *
from utils import units
         
class Cell(dobject):
   """Base class to represent the simulation cell in a periodic system.

   This class has the base attributes required for either flexible or 
   isotropic cell dynamics. Uses an upper triangular matrix to prevent 
   unphysical rotations of the cell (see P. Raiteri, J. Gale and G. Bussi, 
   J. Phys.: Condens. Matter 23 334213 (2011)).

   Attributes:
      h: An array giving the lattice vector matrix.
      p: An array giving the lattice vector momenta.
      m: A float giving the effective cell mass.
      h0: An array giving the reference lattice vector matrix, under no
         external pressure.
      ih: An array giving the inverse of h.
      ih0: An array giving the inverse of h0.
      strain: An array giving the strain tensor.
   """

   def __init__(self, h=None, m=1.0):      
      """Initialises base cell class.

      Args:
         h: Optional array giving the initial lattice vector matrix. The 
            reference cell matrix is set equal to this. Must be an upper
            triangular 3*3 matrix. Defaults to a 3*3 identity matrix.
         m: Optional float giving the effective mass of the cell. Defaults to 1.
      """
      
      if h is None:
         h = np.identity(3, float)      
      dset(self,"h",depend_array(name = 'h', value = h) )
      dset(self,"p",depend_array(name = 'p', value = np.zeros((3,3),float)) )
      dset(self,"m",depend_value(name = "m", value = m) )

      h0 = np.array(h,copy=True)
      dset(self,"h0",depend_array(name = 'h0', value = h0 ) )
            
      dset(self, "ih" , depend_array(name = "ih", value = np.zeros((3,3),float), func=self.get_ih, dependencies=[dget(self,"h")]) )
      dset(self, "ih0" , depend_array(name = "ih0", value = np.zeros((3,3),float), func=self.get_ih0, dependencies=[dget(self,"h0")]) )
      dset(self, "strain", depend_value(name = "strain", func=self.get_strain, dependencies=[dget(self,"h"),dget(self,"ih0")]) )
            
   def get_ih(self):
      """Inverts the lattice vector matrix."""

      return invert_ut3x3(self.h)  
    
   def get_ih0(self):
      """Inverts the reference lattice vector matrix."""

      return invert_ut3x3(self.h0)      

   def get_strain(self):
      """Computes the strain tensor from the unit cell and reference cell."""

      root = np.dot(depstrip(self.h), depstrip(self.ih0))
      eps = np.dot(np.transpose(root), root) - np.identity(3, float)
      eps *= 0.5
      return eps
      
   def apply_pbc(self, atom):
      """Uses the minimum image convention to return a particle to the
         unit cell.

      Args:
         atom: An Atom object.

      Returns:
         An array giving the position of the image that is inside the 
         system box.
      """

      s = np.dot(self.ih,atom.q)
      for i in range(3):
         s[i] = s[i] - round(s[i])
      return np.dot(self.h,s)

   def minimum_distance(self, atom1, atom2):
      """Takes two atoms and tries to find the smallest vector between two 
      images. 

      This is only rigorously accurate in the case of a cubic cell,
      but gives the correct results as long as the cut-off radius is defined
      as smaller than the smallest width between parallel faces.

      Args:
         atom1: An Atom object.
         atom2: An Atom object.

      Returns:
         An array giving the minimum distance between the positions of atoms 
         atom1 and atom2 in the minimum image convention.
      """

      s = np.dot(self.ih,atom1.q-atom2.q)
      for i in range(3):
         s[i] -= round(s[i])
      return np.dot(self.h, s)


class CellFlexi(Cell):
   """Cell object for flexible cell simulations.

   Note that the h matrix must be redefined in the init, so that the upper
   triangular cell matrix can be synchronized with the vector form, and as such
   much of the dependency network must be reworked. The reference lattice vector
   matrix can be defined in this class separately to the system box, so that 
   the initial system box can be in a strained configuration while the 
   reference cell is for the unstrained case, as required.

   Attributes:
      h6: An array giving the 6 non-zero components of the lattice vector 
         matrix in vector form.
      p6: An array giving the 6 non-zero components of the lattice momentum 
         matrix in vector form.
      m6: An array with 6 components all equal to the cell mass. Used when all
         the cell degrees of freedom need to be divided by the cell mass.
      V: A float giving the volume of the system box.
      V0: A float giving the volume of the reference cell.
      kin: The kinetic energy of the cell.
   """

   def __init__(self, h=None, h0=None, m=1.0):    
      """Initialises flexible cell class.

      Args:
         h: Optional array giving the initial lattice vector matrix. The 
            reference cell matrix is set equal to this. Must be an upper
            triangular 3*3 matrix. Defaults to a 3*3 identity matrix.
         h0: Optional array giving the reference lattice vector matrix. 
            Must be an upper triangular 3*3 matrix. Defaults to a h.
         m: Optional float giving the effective mass of the cell. Defaults to 1.
      """

      if h is None:
         h = np.identity(3, float)
      super(CellFlexi,self).__init__(h=h, m=m)

      sync_h=synchronizer()
      sync_p=synchronizer()
      dset(self, "h6", depend_array(name="h6", value=np.zeros(6,float), func={"h":self.htoh6}, synchro=sync_h) )      
      dset(self, "p6", depend_array(name="p6", value=np.zeros(6,float), func={"p":self.ptop6}, synchro=sync_p) )
                
      dset(self, "h", depend_array(name = 'h', value = h, func={"h6":self.h6toh}, synchro=sync_h) )
      dset(self, "p", depend_array(name = 'p', value = np.zeros((3,3),float), func={"p6":self.p6top}, synchro=sync_p) )
      dset(self, "m6", depend_array(name= "m6", value=np.zeros(6,float), func=self.mtom6, dependencies=[dget(self,"m")]) )

      dset(self, "ih", depend_array(name = "ih", value = np.zeros((3,3),float), func=self.get_ih, dependencies=[dget(self,"h")]) )
      dset(self, "strain", depend_value(name = "strain", func=self.get_strain, dependencies=[dget(self,"h"),dget(self,"ih0")]) )
      
      if not h0 is None:
         self.h0 = h0
      
      dset(self, "V", depend_value(name = 'V', func=self.get_volume, dependencies=[dget(self,"h")]) )
      dset(self, "V0", depend_value(name = 'V0', func=self.get_volume0, dependencies=[dget(self,"h0")]) )
      
      dset(self, "kin", depend_value(name = "kin", func=self.get_kin, dependencies=[dget(self,"p6"),dget(self,"m")]) )
      
   def htoh6(self):
      """Transforms the lattice vector matrix from matrix to vector form.

      Returns:
         A vector with 6 components.
      """

      h6 = np.zeros(6, float)
      h = depstrip(self.h)
      h6[0:3] = h[0,0:3]
      h6[3:5] = h[1,1:3]
      h6[5:6] = h[2,2]
      return h6

   def h6toh(self):
      """Transforms the lattice vector matrix from vector to matrix form.

      Returns:
         An upper triangular 3*3 matrix.
      """

      h = np.zeros((3,3), float)
      h6 = depstrip(self.h6)
      h[0,0:3] = h6[0:3]
      h[1,1:3] = h6[3:5]
      h[2,2] = h6[5:6]
      return h

   def ptop6(self):
      """Transforms the lattice momentum matrix from matrix to vector form.

      Returns:
         A vector with 6 components.
      """

      p6 = np.zeros(6, float)
      p = depstrip(self.p)
      p6[0:3] = p[0,0:3]
      p6[3:5] = p[1,1:3]
      p6[5:6] = p[2,2]
      return p6

   def p6top(self):
      """Transforms the lattice momentum matrix from vector to matrix form.

      Returns:
         An upper triangular 3*3 matrix.
      """

      p = np.zeros((3,3), float)
      p6 = depstrip(self.p6)
      p[0,0:3] = p6[0:3]
      p[1,1:3] = p6[3:5]
      p[2,2] = p6[5:6]
      return p
   
   def mtom6(self):
      """Creates a mass vector.

      Makes a 6 component vector with all the components equal to the cell mass.
      Used in algorithms that require that the mass associated with each degree
      of freedom be given, such as in the thermostat.

      Returns:
         A vector with 6 components.
      """

      m6 = np.zeros(6, float)
      m6 = self.m
      return m6
      
   def get_volume(self):
      """Calculates the volume of the system box."""         

      return det_ut3x3(self.h)

   def get_volume0(self):
      """Calculates the volume of the reference box."""

      return det_ut3x3(self.h0)
            
   def get_kin(self):
      """Calculates the kinetic energy of the cell from the cell parameters."""

      p6=depstrip(self.p6)
      return np.dot(p6,p6)/(2.0*self.m)
      
      
class CellRigid(Cell):
   """Cell object for isotropic cell dynamics simulations.

   Note that in this class, as the cell box shape is fixed, the main 
   dynamical variable is the volume, not any of the lattice parameters, so
   the lattice vector matrix is taken to depend upon the volume rather than 
   the other way around.

   Attributes:
      V: A float giving the volume of the system box.
      V0: A float giving the volume of the reference cell.
      P: A float givin the rate of change of volume divided by the cell mass,
         in effect the volume momentum.
      M: An array of one element, containing the mass. Used to access the mass 
         as an array in the thermostating step.
      kin: The kinetic energy of the cell.
   """

   def __init__(self, h=None, m=1.0):    
      """Initialises rigid cell class.

      Args:
         h: Optional array giving the initial lattice vector matrix. The 
            reference cell matrix is set equal to this. Must be an upper
            triangular 3*3 matrix. Defaults to a 3*3 identity matrix.
         m: Optional float giving the effective mass of the cell. Defaults to 1.
      """

      if h is None:
         h=np.identity(3, float)
      super(CellRigid,self).__init__(h, m)

      dset(self, "V0", depend_value(name = 'V0', func=self.get_volume0, dependencies=[dget(self,"h0")]) )
      dset(self, "V", depend_value(name = 'V', value=self.get_volume0()) )

      dset(self, "h", depend_array(name = 'h', value = h, func=self.Vtoh, dependencies=[dget(self,"V"),dget(self,"h0"),dget(self,"V0")]) )
      dset(self, "ih" , depend_array(name = "ih", value = np.zeros((3,3),float), func=self.get_ih, dependencies=[dget(self,"h")]) )

      dset(self, "P", depend_array(name = 'P', value=np.zeros(1,float)) )
      dset(self, "M", depend_array(name="M", value=np.zeros(1,float), func=self.mtoM, dependencies=[dget(self,"m")]) )
      
      dset(self, "p", depend_array(name = 'p', value = np.zeros((3,3),float), func=self.Ptop, dependencies=[dget(self,"P"),dget(self,"V"),dget(self,"h"),dget(self,"ih")]) )

      dset(self, "kin", depend_value(name = "kin", func=self.get_kin, dependencies=[dget(self,"P"),dget(self,"m")]) )
      
   def mtoM(self):
      """Makes the cell mass a vector quantity."""

      return np.identity(1)*self.m

   def Vtoh(self):
      """Gets the lattice vector matrix from the volume."""

      return depstrip(self.h0).copy()*(self.V/self.V0)**(1.0/3.0)

   def Ptop(self):
      """Forms the lattice momentum matrix from the volume momentum."""

      dhdV = np.zeros((3,3))
      dhdV[0,0] = 1.0/(self.h[1,1]*self.h[2,2])
      dhdV[1,1] = 1.0/(self.h[0,0]*self.h[2,2])
      dhdV[2,2] = 1.0/(self.h[1,1]*self.h[0,0])
      return np.dot(dhdV, self.ih)*3*self.V*self.P

   def get_volume0(self):
      """Calculates the volume of the reference box."""

      return det_ut3x3(depstrip(self.h0))
            
   def get_kin(self):
      """Calculates the kinetic energy of the cell from the cell parameters."""

      return self.P[0]**2/(2.0*self.m)

      
class RestartCell(Restart):
   """Cell restart class.

   Handles generating the appropraite cell class from the xml input file,
   and generating the xml checkpoint tags and data from an instance of the 
   object.

   Attributes:
      m: An optional float giving the mass of the cell. Defaults to 0.0.
      h: An optional array giving the system box. Defaults to a 3*3 identity 
         matrix.
      h0: An optional array giving the reference box. Defaults to an array 
         of zeros.
      p: An optional array giving the lattice momentum matrix. Defaults to an 
         array of zeros.
      P: An optional float giving the conjugate momentum to the volume 
         fluctuations. Defaults to 0.0.
      init_temp: An optional float to give the effective temperature that the 
         cell velocities should be initialised to. Defaults to -1.0.
      from_file: An optional string giving the name of a pdb file containing
         the initial cell and atom positions. Defaults to ''.
      flexible: A boolean giving whether the cell will be allowed to change 
         shape. Defaults to False.
   """

   fields={ "m" : (RestartValue, (float, 0.0)), "h" : (RestartArray,(float,np.identity(3))), 
            "h0" : (RestartArray,(float,np.zeros((3,3)))), "p" : (RestartArray,(float,np.zeros((3,3),float))),
            "P" : (RestartValue,(float,0.0)), "init_temp": (RestartValue, (float, -1.0)),
            "from_file": (RestartValue,(str,"")) }
   attribs={ "flexible" : (RestartValue, (bool, False)) }
    
   def __init__(self, cell=None):
      """Initialises RestartCell."""

      super(RestartCell,self).__init__()
      if not cell is None:
         self.store(cell)
      
   def store(self, cell, filename=""):
      """Takes a Cell instance and stores of minimal representation of it.

      Args:
         cell: A cell object.
         filename: An optional float giving a file to read the cell dimensions
            from. Defaults to ''.
      """

      self.from_file.store(filename)
      self.flexible.store(type(cell) is CellFlexi)
      self.m.store(cell.m)
      self.h.store(cell.h)
      self.h0.store(cell.h0)
      if type(cell) is CellFlexi:
         self.p.store(cell.p)
      else:
         self.P.store(cell.P[0])
      self.from_file.store(cell.from_file)
      
   def fetch(self):
      """Creates a cell object.

      Returns: 
         A cell object of the appropriate type and with the appropriate 
         properties given the attributes of the RestartCell object.
      """

      self.check()
      if (self.flexible.fetch()): 
         cell = CellFlexi(h=self.h.fetch(), m=self.m.fetch())
         if det_ut3x3(self.h0.fetch()) == 0.0:
            cell.h0 = cell.h
         cell.p = self.p.fetch()
      else:
         cell = CellRigid(h=self.h.fetch(), m=self.m.fetch())
         cell.P = self.P.fetch()
      cell.from_file = self.from_file.fetch()
      return cell
      
   def check(self):
      """Function that deals with optional arguments.
      
      Deals with the init_temp and from_file arguments, and uses
      them to initialise some of the cell parameters depending on which ones
      have been specified.
      """

      if (self.init_temp.fetch() >= 0):
         pass
      if self.from_file.fetch() != "":
         myatoms, mycell = utils.io.io_pdb.read_pdb(open(self.from_file.fetch(),"r")) 
         self.h.store(mycell.h)
         if (self.h0.fetch() == np.zeros((3,3))).all():
            self.h0.store(mycell.h0)
         if self.m.fetch() == 0.0:
            self.m.store(mycell.m)
         self.from_file.store("")
   
