"""TODO

"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np
import time

from ipi.utils.depend import depend_value, dset, dobject
from ipi.utils.softexit import softexit
from ipi.utils.io import read_file
from ipi.utils.io.inputs.io_xml import xml_parse_file
from ipi.utils.units import unit_to_internal


__all__ = ['Motion', 'ReplayMover']


class Motion(dobject):
    """Base motion calculation class.

    Gives the standard methods and attributes needed in all the
    motion calculation classes.

    Attributes:
        beads: A beads object giving the atoms positions.
        cell: A cell object giving the system box.
        forces: A forces object giving the virial and the forces acting on
            each bead.
        fixcom: A boolean which decides whether the centre of mass
            motion will be constrained or not.
        fixatoms: A list of atoms that should be held fixed to their
            initial positions.

    Depend objects:
        none
    """

    def __init__(self, fixcom=False, fixatoms=None):
        """Initialises Motion object.

        Args:
           fixcom: An optional boolean which decides whether the centre of mass
              motion will be constrained or not. Defaults to False.
           fixatoms: A list of atoms that should be held fixed to their
             initial positions.
        """

        dset(self, "dt", depend_value(name="dt", value=0.0))
        self.fixcom = fixcom
        if fixatoms is None:
            self.fixatoms = np.zeros(0, int)
        else:
            self.fixatoms = fixatoms

    def bind(self, ens, beads, nm, cell, bforce, prng):
        """Binds beads, cell, bforce, and prng to the calculator.

        This takes a beads object, a cell object, a forcefield object and a
        random number generator object and makes them members of the atom motion caclulator.
        It also then creates the objects that will hold the data needed in the
        ensemble algorithms and the dependency network. Note that the conserved
        quantity is defined in the init, but as each ensemble has a different
        conserved quantity the dependencies are defined in bind.

        Args:
            beads: The beads object from whcih the bead positions are taken.
            nm: A normal modes object used to do the normal modes transformation.
            cell: The cell object from which the system box is taken.
            bforce: The forcefield object from which the force and virial are taken.
            bbias: The forcefield object from which the bias forces are obtained
            prng: The random number generator object which controls random number
                generation.
        """

        # store local references to the different bits of the simulation
        self.beads = beads
        self.cell = cell
        self.forces = bforce
        self.prng = prng
        self.nm = nm
        self.ensemble = ens
        self.bias = self.ensemble.bias

    def step(self, step=None):
        """Dummy simulation time step which does nothing."""

        pass


# TODO: Put it in a separate file like the other "movers"?
class ReplayMover(Motion):
    """Calculator object that just loads snapshots from an external file in sequence.

    Has the relevant conserved quantity and normal mode propagator for the
    constant energy ensemble. Note that a temperature of some kind must be
    defined so that the spring potential can be calculated.

    Attributes:
        intraj: The input trajectory file.
        ptime: The time taken in updating the velocities.
        qtime: The time taken in updating the positions.
        ttime: The time taken in applying the thermostat steps.

    Depend objects:
        None really meaningful.
    """

    def __init__(self, fixcom=False, fixatoms=None, intraj=None):
        """Initialises ReplayMover.

        Args:
           dt: The simulation timestep.
           temp: The system temperature.
           fixcom: An optional boolean which decides whether the centre of mass
              motion will be constrained or not. Defaults to False.
           intraj: The input trajectory file.
        """

        super(ReplayMover, self).__init__(fixcom=fixcom, fixatoms=fixatoms)
        if intraj is None:
            raise ValueError("Must provide an initialized InitFile object to read trajectory from")
        self.intraj = intraj
        if intraj.mode == "manual":
            raise ValueError("Replay can only read from PDB or XYZ files -- or a single frame from a CHK file")
        self.rfile = open(self.intraj.value, "r")
        self.rstep = 0

    def step(self, step=None):
        """Does one replay time step."""

        self.ptime = 0.0
        self.ttime = 0.0
        self.qtime = -time.time()

        while True:
            self.rstep += 1
            try:
                if self.intraj.mode == "xyz":
                    for b in self.beads:
                        myatoms = read_file("xyz", self.rfile)
                        myatoms.q *= unit_to_internal("length", self.intraj.units, 1.0)
                        b.q[:] = myatoms.q
                elif self.intraj.mode == "pdb":
                    for b in self.beads:
                        myatoms, mycell = read_file("pdb", self.rfile)
                        myatoms.q *= unit_to_internal("length", self.intraj.units, 1.0)
                        mycell.h *= unit_to_internal("length", self.intraj.units, 1.0)
                        b.q[:] = myatoms.q
                    self.cell.h[:] = mycell.h
                elif self.intraj.mode == "chk" or self.intraj.mode == "checkpoint":

                    # TODO: Adapt the new `Simulation.load_from_xml`?

                    # reads configuration from a checkpoint file
                    xmlchk = xml_parse_file(self.rfile)   # Parses the file.

                    from ipi.inputs.simulation import InputSimulation
                    simchk = InputSimulation()
                    simchk.parse(xmlchk.fields[0][1])
                    mycell = simchk.cell.fetch()
                    mybeads = simchk.beads.fetch()
                    self.cell.h[:] = mycell.h
                    self.beads.q[:] = mybeads.q
                    softexit.trigger(" # Read single checkpoint")
            except EOFError:
                softexit.trigger(" # Finished reading re-run trajectory")
            if (step is None) or (self.rstep > step):
                break

        self.qtime += time.time()
