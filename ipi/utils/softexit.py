"""Utility functions for killing the wrapper softly.

Copyright (C) 2013, Joshua More and Michele Ceriotti

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


Classes:
   Softexit: Concise class to manage cleaning up in case of an emergency exit.
"""

import traceback, sys
from ipi.utils.messages import verbosity, warning

__all__ = ['Softexit', 'softexit']


class Softexit(object):
   """Class to deal with stopping a simulation half way through.

   Holds the functions used to clean up a simulation that has been
   stopped early, either because of a SIGTERM signal or because the
   user has added an EXIT file to the directory in which it is 
   running. This will then properly shut down the socket interface,
   and print out a RESTART file for the appropriate time step.

   Attributes:
      flist: A list of functions used to close down the socket
         interface.
   """

   def __init__(self):
      """Initializes SoftExit."""

      self.flist = []

   def register(self, func):
      """Adds another function to flist.

      Args:
         func: The function to be added to flist.
      """

      self.flist.append(func)

   def trigger(self, message=""):
      """Halts the simulation.

      Prints out a warning message, then runs all the exit functions in flist
      before terminating the simulation.

      Args:
         message: The message to output to standard output.
      """

      if message != "":
         warning("Soft exit has been requested with message: '" + message + "'. Cleaning up.", verbosity.low)
      for f in self.flist:
         f()
      sys.exit()

softexit = Softexit()
