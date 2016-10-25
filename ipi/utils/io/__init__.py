"""Package with functions for reading and writing files.

This module has machinery for abstract I/O handling.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import sys
import os

from ipi.utils.messages import info
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


def read_file(mode, filedesc, **kwargs):
    """Reads one frame from an open `mode`-style file.

    Args:
        filedesc: An open readable file object from a `mode` formatted file.
        All other args are passed directly to the responsible io function.

    Returns:
        A dictionary with 'atoms', 'cell' and 'comment'.
    """

    return importlib.import_module("ipi.utils.io.backends.io_units").\
        process_units(*_get_io_function(mode, "read")(filedesc=filedesc, **kwargs),
                      output=kwargs["output"] if "output" in kwargs.keys() else "objects")


def read_file_name(filename):
    """Read one frame from file, guessing its format from the extension.

    Args:
        filename: Name of input file.

    Returns:
        A dictionary with 'atoms', 'cell' and 'comment'.
    """

    return read_file(os.path.splitext(filename)[1], open(filename))


def iter_file(mode, filedesc):
    """Takes an open `mode`-style file and yields one Atoms object after another.

    Args:
        filedesc: An open readable file object from a `mode` formatted file.

    Returns:
        Generator of frames (dictionary with 'atoms', 'cell' and 'comment') from the trajectory.
    """
    return _get_io_function(mode, "iter")(filedesc=filedesc)


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
            info('Backup performed: {:s} -> {:s}'.format(filename, fn_backup))

    else:
        # There is no need to back up.
        # `open` will sort out whether `mode` is valid.
        pass

    return open(filename, mode, buffering)
