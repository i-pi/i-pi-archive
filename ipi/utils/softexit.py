"""Utility functions for killing the wrapper softly.

Classes:
   Softexit: Concise class to manage cleaning up in case of an emergency exit.
"""

import traceback, sys
from ipi.utils.messages import verbosity, warning

__all__ = ['Softexit', 'softexit']


class Softexit(object):
   def __init__(self):
      self.flist = []

   def register(self, func):
      self.flist.append(func)

   def trigger(self, message=""):

      if message != "":
         warning("Soft exit has been requested with message: '"+message+"'. Cleaning up.", verbosity.low)
      for f in self.flist:
         f()
      sys.exit()

softexit = Softexit()
