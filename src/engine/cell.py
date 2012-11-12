"""Contains the classes which deal with the system box.

Used for implementing the minimum image convention.

Classes:
   Cell: Base cell class with the generic methods and attributes.
"""

__all__ = ['Cell']

import numpy as np
import math
from utils.depend import *
from utils.mathtools import *
from utils import units


class Cell(dobject):
   """Base class to represent the simulation cell in a periodic system.

   This class has the base attributes required for either flexible or
   isotropic cell dynamics. Uses an upper triangular lattice vector matrix to
   represent the cell.

   Depend objects:
      h: An array giving the lattice vector matrix.
   """

   def __init__(self, h=None):
      """Initialises base cell class.

      Args:
         h: Optional array giving the initial lattice vector matrix. The
            reference cell matrix is set equal to this. Must be an upper
            triangular 3*3 matrix. Defaults to a 3*3 identity matrix.
      """

      if h is None:   h = np.identity(3, float)
      dset(self,"h",depend_array(name = 'h', value = h) )

      dset(self,"ih",
         depend_array(name = "ih", value = np.zeros((3,3),float),
            func=self.get_ih, dependencies=[dget(self,"h")]) )
      dset(self,"V",
         depend_value(name = 'V', func=self.get_volume,
            dependencies=[dget(self,"h")]) )

   def get_ih(self):
      """Inverts the lattice vector matrix."""

      return invert_ut3x3(self.h)

   def get_volume(self):
      """Calculates the volume of the system box."""

      return det_ut3x3(self.h)

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

   def array_pbc(self, pos):
      """Uses the minimum image convention to return a list of particles to the
         unit cell.

      Args:
         atom: An Atom object.

      Returns:
         An array giving the position of the image that is inside the
         system box.
      """

      s = pos.copy()
      s.shape=(len(pos)/3,3)

      s = np.dot(self.ih,s.T)
      s = s - np.round(s)

      s=np.dot(self.h,s).T

      pos=s.reshape((len(s)*3))


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
