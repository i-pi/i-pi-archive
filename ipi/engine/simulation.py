"""A class that runs the simulation and outputs results.

The root class for the whole simulation. Contains references to all the top
level objects used in the simulation, and controls all the steps that are
not inherently system dependent, like the running of each time step,
choosing which properties to initialise, and which properties to output.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import os, threading
import time
from copy import deepcopy

from ipi.utils.depend import depend_value, dobject, dset
from ipi.utils.io.inputs.io_xml import xml_parse_file
from ipi.utils.messages import verbosity, info, warning, banner
from ipi.utils.softexit import softexit
import ipi.engine.outputs as eoutputs
import ipi.inputs.simulation as isimulation


__all__ = ['Simulation']


class Simulation(dobject):
    """Main simulation object.

    Contains all the references and the main dynamics loop. Also handles the
    initialisation and output.

    Attributes:
        prng: A random number generator object.
        mode: A string specifying what kind of simulation will be run.
        fflist: A list of forcefield objects that can be called to compute energy and forces
        syslist: A list of physical systems
        prng: A random number generator.
        tsteps: The total number of steps.
        ttime: The wall clock time (in seconds).
        outtemplate: A template output object to be used to generate the outputs
            list. This will be used for each of the systems.
        outputs: A list of output objects that should be printed during the run
        paratemp: A helper object for parallel tempering simulations
        chk: A checkpoint object which is kept up-to-date in case of emergency exit
        rollback: If set to true, the state of the simulation at the start
            of the step will be output to a restart file rather than
            the current state of the simulation. This is because we cannot
            restart from half way through a step, only from the beginning of a
            step, so this is necessary for the trajectory to be continuous.

    Depend objects:
        step: The current simulation step.
    """

    @staticmethod
    def load_from_xml(fn_input, custom_verbosity=None, request_banner=False):
        """Load an XML input file and return a `Simulation` object.

        Arguments:
            fn_input (str): Name of the input file.
            custom_verbosity (str): If not `None`, overrides the verbosity
                specified by the input file.
            request_banner (bool): Whether to print the i-PI banner,
                if verbosity is higher than 'quiet'.
        """

        # parse the file
        xmlrestart = xml_parse_file(open(fn_input))

        # prepare the simulation input object
        input_simulation = isimulation.InputSimulation()

        # check the input and partition it appropriately
        input_simulation.parse(xmlrestart.fields[0][1])

        # override verbosity if requested
        if custom_verbosity is not None:
            input_simulation.verbosity.value = custom_verbosity

        # print banner if not suppressed and simulation verbose enough
        if request_banner and input_simulation.verbosity.fetch() != 'quiet':
            banner()

        # create the simulation object
        simulation = input_simulation.fetch()

        # pipe between the components of the simulation
        simulation.bind()

        # echo the input file if verbose enough
        if verbosity.level > 0:
            print " # i-PI loaded input file: ", fn_input
        if verbosity.level > 1:
            print " --- begin input file content ---"
            ifile = open(fn_input, "r")
            for line in ifile.readlines():
                print line,
            print " ---  end input file content  ---"
            ifile.close()

        return simulation

    def __init__(self, mode, syslist, fflist, outputs, prng, smotion=None, step=0, tsteps=1000, ttime=0, threads=False):
        """Initialises Simulation class.

        Args:
            mode: What kind of simulation is this
            syslist: A list of system objects
            fflist: A list of forcefield objects
            prng: A random number object.
            smotion: A "super-motion" class specifying what to do with different system replicas
            outputs: A list of output objects.
            step: An optional integer giving the current simulation time step.
                Defaults to 0.
            tsteps: An optional integer giving the total number of steps. Defaults
                to 1000.
            ttime: The simulation running time. Used on restart, to keep a
                cumulative total.
        """

        info(" # Initializing simulation object ", verbosity.low)
        self.prng = prng
        self.mode = mode
        self.threading = threads

        self.syslist = syslist
        for s in syslist:
            s.prng = self.prng    # bind the system's prng to self prng
            s.init.init_stage1(s)

        #! TODO - does this have any meaning now that we introduce the smotion class?
        if self.mode == "md" and len(syslist) > 1:
            warning("Multiple systems will evolve independently in a '" + self.mode + "' simulation.")

        self.fflist = {}
        for f in fflist:
            self.fflist[f.name] = f

        self.outtemplate = outputs

        dset(self, "step", depend_value(name="step", value=step))
        self.tsteps = tsteps
        self.ttime = ttime
        self.smotion = smotion

        self.chk = None
        self.rollback = True

    def bind(self):
        """Calls the bind routines for all the objects in the simulation."""

        if self.tsteps <= self.step:
            raise ValueError("Simulation has already run for total_steps, will not even start. "
                             "Modify total_steps or step counter to continue.")

        for s in self.syslist:
            # binds important computation engines
            s.bind(self)

        self.outputs = []
        for o in self.outtemplate:
            if type(o) is eoutputs.CheckpointOutput:    # checkpoints are output per simulation
                o.bind(self)
                self.outputs.append(o)
            else:   # properties and trajectories are output per system
                isys = 0
                for s in self.syslist:   # create multiple copies
                    no = deepcopy(o)
                    if s.prefix != "":
                        no.filename = s.prefix + "_" + no.filename
                    no.bind(s)
                    self.outputs.append(no)
                    isys += 1

        self.chk = eoutputs.CheckpointOutput("RESTART", 1, True, 0)
        self.chk.bind(self)

        if not self.smotion is None:
            self.smotion.bind(self.syslist, self.prng)

    def softexit(self):
        """Deals with a soft exit request.

        Tries to ensure that a consistent restart checkpoint is
        written out.
        """

        if self.step < self.tsteps:
            self.step += 1
        if not self.rollback:
            info("SOFTEXIT: Saving the latest status at the end of the step")
            self.chk.store()

        self.chk.write(store=False)

    def run(self):
        """Runs the simulation.

        Does all the simulation steps, and outputs data to the appropriate files
        when necessary. Also deals with starting and cleaning up the threads used
        in the communication between the driver and the PIMD code.
        """

        # registers the softexit routine
        softexit.register_function(self.softexit)
        softexit.start(self.ttime)

        for k, f in self.fflist.iteritems():
            f.run()


        # prints inital configuration -- only if we are not restarting
        if self.step == 0:
            self.step = -1
            # must use multi-threading to avoid blocking in multi-system runs with WTE
            if self.threading:
                stepthreads = []
                for o in self.outputs:
                    st = threading.Thread(target=o.write, name=o.filename)
                    st.daemon = True
                    st.start()
                    stepthreads.append(st)

                for st in stepthreads:
                    while st.isAlive():
                        # This is necessary as join() without timeout prevents main from receiving signals.
                        st.join(2.0)
            else:
                for o in self.outputs:
                    o.write()  # threaded output seems to cause random hang-ups. should make things properly thread-safe

            self.step = 0

        steptime = 0.0
        simtime = time.time()

        cstep = 0
        #tptime = 0.0
        #tqtime = 0.0
        #tttime = 0.0
        ttot = 0.0
        # main MD loop
        for self.step in xrange(self.step, self.tsteps):
            # stores the state before doing a step.
            # this is a bit time-consuming but makes sure that we can honor soft
            # exit requests without screwing the trajectory

            steptime = -time.time()
            if softexit.triggered:
                break

            self.chk.store()

            if self.threading:
                stepthreads = []
                # steps through all the systems
                for s in self.syslist:
                    # creates separate threads for the different systems
                    st = threading.Thread(target=s.motion.step, name=s.prefix, kwargs={"step":self.step})
                    st.daemon = True
                    st.start()
                    stepthreads.append(st)

                for st in stepthreads:
                    while st.isAlive():
                        # This is necessary as join() without timeout prevents main from receiving signals.
                        st.join(2.0)
            else:
                for s in self.syslist:
                    s.motion.step(step=self.step)
                    
            if softexit.triggered:
                # Don't continue if we are about to exit.
                break
            
            # does the "super motion" step
            if self.smotion is not None:
                # TODO: We need a file where we store the exchanges
                self.smotion.step(self.step)

            if softexit.triggered:
                # Don't write if we are about to exit.
                break
            
            if self.threading:
                stepthreads = []
                for o in self.outputs:
                    st = threading.Thread(target=o.write, name=o.filename)
                    st.daemon = True
                    st.start()
                    stepthreads.append(st)

                for st in stepthreads:
                    while st.isAlive():
                        # This is necessary as join() without timeout prevents main from receiving signals.
                        st.join(2.0)
            else:
                for o in self.outputs:
                    o.write()

            steptime += time.time()
            ttot += steptime
            cstep += 1

            if (verbosity.high or (verbosity.medium and self.step % 100 == 0) or (verbosity.low and self.step % 1000 == 0)):
                info(" # Average timings at MD step % 7d. t/step: %10.5e" % (self.step, ttot / cstep))
                cstep = 0
                ttot = 0.0
                #info(" # MD diagnostics: V: %10.5e    Kcv: %10.5e   Ecns: %10.5e" %
                #     (self.properties["potential"], self.properties["kinetic_cv"], self.properties["conserved"] ) )

            if os.path.exists("EXIT"):
                info(" # EXIT file detected! Bye bye!", verbosity.low)
                break

            if (self.ttime > 0) and (time.time() - simtime > self.ttime):
                info(" # Wall clock time expired! Bye bye!", verbosity.low)
                break

        self.rollback = False
