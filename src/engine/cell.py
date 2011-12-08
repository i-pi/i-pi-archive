import numpy as np
import math
from utils.depend import *
from utils.mathtools import det_ut3x3, invert_ut3x3
from utils import units
import pdb
class Cell(dobject):
   """Represents the simulation cell in a periodic system
      Contains: h = lattice vector matrix, p = lattice momentum matrix (NST),
      pc = lattice scalar momentum (NPT), w = barostat mass, pext = external 
      pressure tensor, V = volume, ih = inverse lattice matrix, 
      (h0, ih0, V0) = as above, but for reference cell (pext = 0), 
      kin = kinetic energy, pot = strain potential, 
      strain = strain tensor, piext = external stress tensor
      Initialised by: cell = Cell(h, w, h0, pext)
      h is the lattice vector matrix
      w = barostat mass, default = 1.0
      h0 = reference cell, default = h
      pext = external pressure tensor, default = 0"""

   def __init__(self, h = numpy.identity(3, float)):      
      #un-dependent properties
      dset(self,"h",depend_array(name = 'h', value = h) )
      dset(self,"p",depend_array(name = 'p', value = numpy.zeros((3,3),float)) )
      dset(self,"V",depend_value(name = 'V', deps=depend_func(func=self.get_volume, dependencies=[depget(self,"h")])) )
      dset(self,"m",depend_value(name = "m", value = 1.0) )
      
      dset(self,"kin", depend_value(name = "kin", deps=depend_func(func=self.get_kin, dependencies=[depget(self,"p"),depget(self,"m")])) )
      dset(self,"ih" , depend_array(name = "ih", value = numpy.zeros((3,3),float), deps=depend_func(func=self.get_ih, dependencies=[depget(self,"h")])) )
      
      
   def get_volume(self):
      """Calculates the volume of the unit cell, assuming an upper-triangular
         lattice vector matrix"""         
      return det_ut3x3(self.h)

   def get_kin(self):
      """Calculates the kinetic energy of the cell from the cell parameters"""
      p=np.array(self.p)   # strips the deps object from p
      ke = 0.0
      for i in range(3):
         for j in range(i,3):
            ke += p[i, j]**2            
      ke /= 2.0*self.m
      return ke
      
   def get_ih(self):
      """Inverts a 3*3 (upper-triangular) cell matrix"""
      return invert_ut3x3(self.h)      

   def apply_pbc(self, atom):
      """Uses the minimum image convention to return a particle to the
         unit cell"""

      s=numpy.dot(self.ih,atom.q)
      for i in range(3):
         s[i] = s[i] - round(s[i])
      return numpy.dot(self.h,s)

   def minimum_distance(self, atom1, atom2):
      """Takes two atoms and tries to find the smallest vector between two 
         images. This is only rigorously accurate in the case of a cubic cell,
         but gives the correct results as long as the cut-off radius is defined
         as smaller than the smallest width between parallel faces."""

      s=numpy.dot(self.ih,atom1.q-atom2.q)
      
      for i in range(3):
         s[i] -= round(s[i])
         
      return numpy.dot(self.h, s)

