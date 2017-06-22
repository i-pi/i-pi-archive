"""Package with functions for reading and writing files.

This module has machinery for abstract I/O handling.

The idea is that the unit conversion is done here. The default is to guess
units from the file, but it can be overridden. 
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import sys
import os

from ipi.utils.messages import info, verbosity
from ipi.utils.units import unit_to_user
from ipi.external import importlib
from ipi.utils.decorators import cached

__all__ = ["iter_file", "print_file_path", "print_file", "read_file"]


mode_map = {
    "bin":  "binary",
}


io_map = {
    "print_path": "print_%s_path",
    "print":      "print_%s",
    "read":       "read_%s",
    "iter":       "iter_%s",
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

    try:
        mode = mode[mode.find(".")+1:]
        if mode in mode_map:
            mode = mode_map[mode]
        module = importlib.import_module("ipi.utils.io.backends.io_%s" % mode)
    except ImportError:
        print "Error: mode %s is not supported." % mode
        sys.exit()

    try:
        func = getattr(module, io_map[io] % mode)
    except KeyError:
        print "Error: io %s is not supported with mode %s." % (io, mode)
        sys.exit()

    return func

# VENKAT TODO: must overhaul also this function
def print_file_path(mode, beads, cell, filedesc=sys.stdout, title="", key="", dimension="length", units="automatic", cell_units="automatic"):
    """Prints all the bead configurations, into a `mode` formatted file.

    Prints all the replicas for each time step separately, rather than all at
    once.

    Args:
        beads: A beads object giving the bead positions.
        cell: A cell object giving the system box.
        filedesc: An open writable file object. Defaults to standard output.
    """
    
    return _get_io_function(mode, "print_path")(beads=beads, cell=cell, filedesc=filedesc)

# VENKAT TODO : also get the print_file functions work with just arrays, so we have a "fast write" mode that sidesteps any parsing or conversion, similar to what I'm doing for readfile and readfile_raw

def print_file_raw(mode, atoms, cell, filedesc=sys.stdout, title="", key="", dimension="length", units="automatic", cell_units="automatic"):
    """Prints the centroid configurations, into a `mode` formatted file.

    Args:
        atoms: An atoms object giving the centroid positions.
        cell: A cell object giving the system box.
        filedesc: An open writable file object. Defaults to standard output.
        title: This gives a string to be appended to the comment line.
    """
   
    return _get_io_function(mode, "print")(atoms=atoms, cell=cell, filedesc=filedesc, title=title, cell_conv=cell_conv, atoms_conv=atoms_conv)

def print_file(mode, atoms, cell, filedesc=sys.stdout, title="", key="", dimension="length", units="automatic", cell_units="automatic"):
    """Prints the centroid configurations, into a `mode` formatted file.

    Args:
        atoms: An atoms object giving the centroid positions.
        cell: A cell object giving the system box.
        filedesc: An open writable file object. Defaults to standard output.
        title: This gives a string to be appended to the comment line.
    """
 
    if mode == "pdb":   # special case for PDB
        if dimension != "length":
            raise ValueError("PDB Standard is only designed for atomic positions")
        if units == "automatic": units="angstrom"
        if cell_units == "automatic": cell_units="angstrom"
    # in general, "automatic" units are actually "atomic_units"
    else:
        if units == "automatic": units="atomic_unit"
        if cell_units == "automatic": cell_units="atomic_unit"
 
    cell_conv = unit_to_user("length", cell_units, 1.0)
    atoms_conv = unit_to_user(dimension, units, 1.0)
 
    title = title + ("%s{%s}  cell{%s}" % (key, units, cell_units))

    print_file_raw(mode=mode, atoms=atoms, cell=cell, filedesc=filedesc, title=title, cell_conv=cell_conv, atoms_conv=atoms_conv)

def read_file_raw(mode, filedesc):    

    reader = _get_io_function(mode, "read")
        
    comment, cell, atoms, names, masses = reader(filedesc=filedesc)
     
    return {
          "comment" : comment, 
          "data": atoms,
          "masses": masses,
          "names": names,
          "natoms": len(names),
          "cell": cell
        }
 

def read_file(mode, filedesc, dimension="", units="automatic", cell_units="automatic"):
    """Reads one frame from an open `mode`-style file.

    Args:
        filedesc: An open readable file object from a `mode` formatted file.
        All other args are passed directly to the responsible io function.

    Returns:
        A dictionary as returned by `process_units`.
    """

    print "reading file", dimension, units, cell_units
    raw_read = read_file_raw(mode=mode, filedesc=filedesc)

    # late import is needed to break an import cycle
    from .io_units import process_units
    
    return process_units(dimension=dimension, units=units, cell_units=cell_units, mode=mode, **raw_read)


def read_file_name(filename):
    """Read one frame from file, guessing its format from the extension.

    Args:
        filename: Name of input file.

    Returns:
        A dictionary with 'atoms', 'cell' and 'comment'.
    """

    return read_file(os.path.splitext(filename)[1], open(filename))

# VENKAT TODO also to the same on iter_file, and get rid of this useless kwargs thing
def iter_file(mode, filedesc, output="objects", **kwargs):
    """Takes an open `mode`-style file and yields one Atoms object after another.

    Args:
        filedesc: An open readable file object from a `mode` formatted file.

    Returns:
        Generator of frames dictionaries, as returned by `process_units`.
    """

    # late import is needed to break an import cycle
    from .io_units import process_units

    reader = _get_io_function(mode, "read")

    try:
        while True:
            yield process_units(*reader(filedesc=filedesc, **kwargs), output=output)
    except EOFError:
        pass


def iter_file_name(filename):
    """Open a trajectory file, guessing its format from the extension.

    Args:
        filename: Filename of a trajectory file.

    Returns:
        Generator of frames (dictionary with 'atoms', 'cell' and 'comment') from the trajectory in `filename`.
    """

    return iter_file(os.path.splitext(filename)[1], open(filename))


def open_backup(filename, mode='r', buffering=-1):
    """A wrapper around `open` which saves backup files.

    If the file is opened in write mode and already exists, it is first
    backed up under a new file name, keeping all previous backups. Then,
    a new file is opened for writing.

    For reference: https://docs.python.org/2/library/functions.html#open

    Args:
        The same as for `open`.

    Returns:
        An open file as returned by `open`.
    """

    if mode.startswith('w'):

        # If writing, make sure nothing is overwritten.

        i = 0
        fn_backup = filename
        while os.path.isfile(fn_backup):
            fn_backup = '#' + filename + '#%i#' % i
            i += 1

        if fn_backup != filename:
            os.rename(filename, fn_backup)
            info('Backup performed: {:s} -> {:s}'.format(filename, fn_backup), verbosity.low)

    else:
        # There is no need to back up.
        # `open` will sort out whether `mode` is valid.
        pass

    return open(filename, mode, buffering)
