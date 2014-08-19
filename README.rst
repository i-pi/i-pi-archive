====
i-PI
====

A Python interface for ab initio path integral molecular dynamics simulations.
i-PI is composed of a Python server (i-pi itself, that does not need to be
compiled but only requires a relatively recent version of Python and Numpy)
that propagates the (path integral) dynamics of the nuclei, and of an external
code that acts as client and computes the electronic energy and forces.

This is typically a patched version of an electronic structure code, but a
simple self-contained Fortran driver that implements Lennard-Jones and
Silvera-Goldman potentials is included for test purposes.


Quick Installation and Test
===========================

Follow these instructions to test i-PI. It is assumed that i-PI will be run
from a Linux environment, with a recent version of Python, Numpy and gfortran,
and that the terminal is initially in the i-pi package directory (the directory
containing this file).

Compile the driver code
-----------------------

::

  $ cd driver
  $ make
  $ cd ..

Run one of the examples
-----------------------

This will first start the wrapper in the background, redirecting the output to
a log file, and then run a couple of instances of the driver code. The progress
of the wrapper is followed by monitoring the log file with the `tail` Linux
command.

::

  $ cd examples/tutorial/tutorial-1/
  $ ../../../i-pi tutorial-1.xml > log &
  $ ../../../drivers/driver.x -h localhost -p 31415 -m sg -o 15 &
  $ ../../../drivers/driver.x -h localhost -p 31415 -m sg -o 15 &
  $ tail -f log

The monitoring can be interrupted with CTRL+C when the run has finished (5000 steps).

Run the automatic test suite
----------------------------

The automatic test suite can be run with the python package `nose` from the
root directory of the i-pi package.

::

  $ nosetests -v
