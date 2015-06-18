"""Contains the functions used to print the trajectories and read input
configurations with arbitary formatting.

Copyright (C) 2014, i-PI Developer Team

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

Functions:
   print_file_path: Prints all the bead configurations.
   print_file: Prints the centroid configurations.
   read_file: Reads the cell parameters and atom configurations from a xyz file.
   iter_file: Iterate over cell parameters and atom configurations from a xyz file.
"""
import sys

from ipi.utils.decorators import cached

# For backwards compatibility import old known functions from backend and inputs.
from ipi.utils.io.backends import io_pdb, io_xyz, io_binary
from ipi.utils.io.inputs import io_xml

__all__ = [ "io_xml", "io_pdb" , "io_xyz", "io_binary" ]

mode_map = {
         "bin"   :  "binary",
      }

io_map = {
         "print_path"   :  "print_%s_path",
         "print"        :  "print_%s",
         "read"         :  "read_%s",
         "iter"         :  "iter_%s",
      }

@cached
def _get_io_function(mode, io):
   """Returns io function with specified mode.

   This will be determined on the fly based on file name extension and
   available ipi/utils/io/backends/io_*.py backends.

   Args:
      mode: Which format has the file? e.g. "pdb", "xml" or "xyz"
      io: One of "print_path", "print", "read" or "iter"
   """
   import importlib

   try:
      mode = mode[mode.find(".")+1:]
      if mode in mode_map:
         mode = mode_map[mode]
      module = importlib.import_module("ipi.utils.io.backends.io_%s"%mode)
   except ImportError:
      print "Error: mode %s is not supported."% mode
      sys.exit()

   try:
      func = getattr(module, io_map[io]%mode)
   except KeyError:
      print "Error: io %s is not supported with mode %s."% (io,mode)
      sys.exit()

   return func

def print_file_path(mode, beads, cell, filedesc=sys.stdout):
   """Prints all the bead configurations, into a `mode` formatted file.

   Prints all the replicas for each time step separately, rather than all at
   once.

   Args:
      beads: A beads object giving the bead positions.
      cell: A cell object giving the system box.
      filedesc: An open writable file object. Defaults to standard output.
   """
   return _get_io_function(mode, "print_path")(beads=beads, cell=cell, filedesc=filedesc)

def print_file(mode, atoms, cell, filedesc=sys.stdout, title=""):
   """Prints the centroid configurations, into a `mode` formatted file.

   Args:
      atoms: An atoms object giving the centroid positions.
      cell: A cell object giving the system box.
      filedesc: An open writable file object. Defaults to standard output.
      title: This gives a string to be appended to the comment line.
   """
   return _get_io_function(mode, "print")(atoms=atoms, cell=cell, filedesc=filedesc, title=title)

def read_file(mode, filedesc):
   """Takes a `mode`-style file and creates an Atoms object.

   Args:
      filedesc: An open readable file object from a `mode` formatted file.

   Returns:
      An Atoms object with the appropriate atom labels, masses and positions.
   """
   return _get_io_function(mode, "read")(filedesc=filedesc)

def iter_file(mode, filedesc):
   """Takes a `mode`-style file and yields one Atoms object after another.

   Args:
      filedesc: An open readable file object from a `mode` formatted file.

   Returns:
      Generator over the xyz trajectory, that yields
      Atoms objects with the appropriate atom labels, masses and positions.
   """
   return _get_io_function(mode, "iter")(filedesc=filedesc)
