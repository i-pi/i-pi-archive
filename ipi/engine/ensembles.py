"""Contains the classes that deal with the different dynamics required in
different types of ensembles.

Holds the algorithms required for normal mode propagators, and the objects to
do the constant temperature and pressure algorithms. Also calculates the
appropriate conserved energy quantity for the ensemble of choice.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import time

import numpy as np

from ipi.utils.depend import *
from ipi.utils import units
from ipi.utils.softexit import softexit
from ipi.utils.io.backends.io_xyz import read_xyz
from ipi.utils.io.backends.io_pdb import read_pdb
from ipi.utils.io.inputs.io_xml import xml_parse_file
from ipi.utils.units import Constants, unit_to_internal
from ipi.inputs.thermostats import InputThermo
from ipi.inputs.barostats import InputBaro
from ipi.engine.thermostats import *
from ipi.engine.barostats import *


__all__ = ['Ensemble']


class Ensemble(dobject):
   """Base ensemble class.

   Defines the thermodynamic state of the system.   

   Depend objects:
      temp: The system's temperature.
      pext: The systems's pressure
      stressext: The system's stress tensor
   """

   def __init__(self, eens = 0.0, temp = None, pext = None, stressext = None):
      """Initialises Ensemble.

      Args:
         dt: The timestep of the simulation algorithms.
         temp: The temperature.
         fixcom: An optional boolean which decides whether the centre of mass
            motion will be constrained or not. Defaults to False.
      """
      
      dset(self, "temp",  depend_value(name='temp'))
      if not temp is None:
         self.temp = temp
      else: self.temp =0.0
      
      dset(self,"stressext",depend_array(name='stressext',value=np.zeros((3,3),float)))
      if not stressext is None:
         self.stressext = stressext
      else: self.stressext = 0.0

      dset(self,"pext",depend_value(name='pext'))
      if not pext is None:
         self.pext = pext
      else: self.pext = 0.0
      
      dset(self, "eens",  depend_value(name='eens'))
      if not eens is None:
         self.eens = eens
      else: self.eens =0.0      

