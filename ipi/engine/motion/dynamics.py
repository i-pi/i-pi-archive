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

from ipi.engine.motion import Motion
from ipi.utils.depend import depstrip, depend_value, dget, dset, dobject, deppipe
from ipi.engine.thermostats import Thermostat
from ipi.engine.barostats import Barostat

#__all__ = ['Dynamics', 'NVEIntegrator', 'NVTIntegrator', 'NPTIntegrator', 'NSTIntegrator', 'SCIntegrator`']

class Dynamics(Motion):
    """self (path integral) molecular dynamics class.

    Gives the standard methods and attributes needed in all the
    dynamics classes.

    Attributes:
        beads: A beads object giving the atoms positions.
        cell: A cell object giving the system box.
        forces: A forces object giving the virial and the forces acting on
            each bead.
        prng: A random number generator object.
        nm: An object which does the normal modes transformation.

    Depend objects:
        econs: The conserved energy quantity appropriate to the given
            ensemble. Depends on the various energy terms which make it up,
            which are different depending on the ensemble.he
        temp: The system temperature.
        dt: The timestep for the algorithms.
        ntemp: The simulation temperature. Will be nbeads times higher than
            the system temperature as PIMD calculations are done at this
            effective classical temperature.
    """

    def __init__(self, timestep, mode="nve", splitting="obabo", thermostat=None, barostat=None, fixcom=False, fixatoms=None, nmts=None):
        """Initialises a "dynamics" motion object.

        Args:
            dt: The timestep of the simulation algorithms.
            fixcom: An optional boolean which decides whether the centre of mass
                motion will be constrained or not. Defaults to False.
        """

        super(Dynamics, self).__init__(fixcom=fixcom, fixatoms=fixatoms)

        # initialize time step. this is the master time step that covers a full time step
        dset(self, "dt", depend_value(name='dt', value=timestep))
                
        if thermostat is None:
            self.thermostat = Thermostat()
        else:
            self.thermostat = thermostat

        if barostat is None:
            self.barostat = Barostat()
        else:
            self.barostat = barostat

        # multiple time stepping array. 
        if nmts is None or len(nmts)==0:
           self.nmts = np.asarray([1],int)      
        else:
           self.nmts=np.asarray(nmts)
        
        self.enstype = mode
        if self.enstype == "nve":
            self.integrator = NVEIntegrator()
        elif self.enstype == "nvt":
            self.integrator = NVTIntegrator()
        elif self.enstype == "npt":
            self.integrator = NPTIntegrator()
        elif self.enstype == "nst":
            self.integrator = NSTIntegrator()
        elif self.enstype == "sc":
            self.integrator = SCIntegrator()        
        else:
            self.integrator = DummyIntegrator()

        # splitting mode for the integrators
        self.splitting = splitting        
        
        # constraints
        self.fixcom = fixcom
        if fixatoms is None:
            self.fixatoms = np.zeros(0, int)
        else:
            self.fixatoms = fixatoms

    def bind(self, ens, beads, nm, cell, bforce, prng):
        """Binds ensemble beads, cell, bforce, and prng to the dynamics.

        This takes a beads object, a cell object, a forcefield object and a
        random number generator object and makes them members of the ensemble.
        It also then creates the objects that will hold the data needed in the
        ensemble algorithms and the dependency network. Note that the conserved
        quantity is defined in the init, but as each ensemble has a different
        conserved quantity the dependencies are defined in bind.

        Args:
            beads: The beads object from whcih the bead positions are taken.
            nm: A normal modes object used to do the normal modes transformation.
            cell: The cell object from which the system box is taken.
            bforce: The forcefield object from which the force and virial are
                taken.
            prng: The random number generator object which controls random number
                generation.
        """

        super(Dynamics, self).bind(ens, beads, nm, cell, bforce, prng)

        # n times the temperature (for the path integral partition function)
        dset(self, "ntemp", depend_value(name='ntemp', func=self.get_ntemp,
             dependencies=[dget(self.ensemble, "temp")]))

        # fixed degrees of freedom count
        fixdof = len(self.fixatoms) * 3 * self.beads.nbeads
        if self.fixcom:
            fixdof += 3

        # first makes sure that the thermostat has the correct temperature, then proceed with binding it.
        deppipe(self, "ntemp", self.thermostat, "temp")

        # depending on the kind, the thermostat might work in the normal mode or the bead representation.
        self.thermostat.bind(beads=self.beads, nm=self.nm, prng=prng, fixdof=fixdof)

        # binds the barostat
        deppipe(self, "ntemp", self.barostat, "temp")        
        deppipe(self.ensemble, "pext", self.barostat, "pext")
        deppipe(self.ensemble, "stressext", self.barostat, "stressext")
        
        self.barostat.bind(beads, nm, cell, bforce, prng=prng, fixdof=fixdof)

        # now we need to define timesteps for the different propagators. 
        # this depends on the splittings we'll do. thanks to depend machinery
        # everything will be propagated down to where it's needed        
        dset(self, "halfdt", depend_value(name='dt', func=(lambda : self.dt*0.5), dependencies=[dget(self,"dt")]) )

        # O=stochastic B=momenta A=positions 
        if self.splitting == "obabo" :            
            deppipe(self, "halfdt", self.thermostat, "dt")
            deppipe(self, "halfdt", self.barostat, "dt")
                
            # the free ring polymer propagator is called in the inner loop, so propagation time should be redefined accordingly. 
            self.inmts = 1
            for mk in self.nmts: self.inmts*=mk            
            dset(self,"deltat", depend_value(name="deltat", func=(lambda : self.dt/self.inmts) , dependencies=[dget(self,"dt")]) )
            deppipe(self,"deltat", self.nm, "dt")
        elif self.splitting == "baoab" :                        
            deppipe(self, "halfdt", self.barostat, "dt")
                
            # the free ring polymer propagator is called in the inner loop, so propagation time should be redefined accordingly. 
            self.inmts = 1
            for mk in self.nmts: self.inmts*=mk
            dset(self,"deltat", depend_value(name="deltat", func=(lambda : self.dt/self.inmts) , dependencies=[dget(self,"dt")]) )
            dset(self,"halfdeltat", depend_value(name="halfdeltat", func=(lambda : 0.5*self.dt/self.inmts) , dependencies=[dget(self,"dt")]) )                
            
            deppipe(self, "deltat", self.thermostat, "dt")
            if self.enstype == "nve":
                deppipe(self,"deltat", self.nm, "dt")
            else:
                deppipe(self,"halfdeltat", self.nm, "dt")
        elif self.splitting == "aboba" :            
            deppipe(self, "dt", self.thermostat, "dt")
            deppipe(self, "halfdt", self.barostat, "dt")
                
            # the free ring polymer propagator is called in the inner loop, so propagation time should be redefined accordingly. 
            if len(self.nmts)>1 or self.nmts[0]>1:
                raise ValueError("General MTS is not implemented with ABOBA integrator")
            deppipe(self,"halfdt", self.nm, "dt")
        
        self.ensemble.add_econs(dget(self.thermostat, "ethermo"))
        self.ensemble.add_econs(dget(self.barostat, "ebaro"))

        #!TODO THOROUGH CLEAN-UP AND CHECK
        #if self.enstype in ["nvt", "npt", "nst"]:
        if self.enstype == "nvt" or self.enstype == "npt" or self.enstype == "nst":
            if self.ensemble.temp < 0:
                raise ValueError("Negative or unspecified temperature for a constant-T integrator")
            if self.enstype == "npt":
                if type(self.barostat) is Barostat:
                    raise ValueError("The barostat and its mode have to be specified for constant-p integrators")
                if self.ensemble.pext < 0:
                    raise ValueError("Negative or unspecified pressure for a constant-p integrator")
            elif self.enstype == "nst":
                if len(self.nmts)>1 or self.nmts[0]>1:
                    raise ValueError("MTS is not implemented yet for Parrinello-Rahman integrator")
                if np.trace(self.ensemble.stressext) < 0:
                    raise ValueError("Negative or unspecified stress for a constant-s integrator")

		# Binds integrators
        self.integrator.bind(self)
        self.integrator.pconstraints()


    def get_ntemp(self):
        """Returns the PI simulation temperature (P times the physical T)."""

        return self.ensemble.temp * self.beads.nbeads

    def step(self, step=None):
        self.integrator.step(step)        


class DummyIntegrator(dobject):
    """ No-op integrator for (PI)MD """

    def __init__(self):

        pass

    def bind(self, motion):
        """ Reference all the variables for simpler access."""

        self.beads = motion.beads
        self.bias = motion.bias
        self.ensemble = motion.ensemble
        self.forces = motion.forces
        self.prng = motion.prng
        self.nm = motion.nm
        self.thermostat = motion.thermostat
        self.barostat = motion.barostat
        self.fixcom = motion.fixcom
        self.fixatoms = motion.fixatoms
        self.splitting = motion.splitting
        dset(self, "dt", dget(motion, "dt"))
        dset(self, "halfdt", dget(motion, "halfdt"))
        self.nmts = motion.nmts        
        #mts on sc force in suzuki-chin
        if motion.enstype == "sc":
            if(motion.nmts.size > 1):
                raise ValueError("MTS for SC is not implemented yet....")
            else:
                # coefficients to get the (baseline) trotter to sc conversion
                self.coeffsc = np.ones((self.beads.nbeads,3*self.beads.natoms), float)
                self.coeffsc[::2] /= -3.
                self.coeffsc[1::2] /= 3.                
                self.nmts=motion.nmts[-1]

    def pstep(self):
        """Dummy momenta propagator which does nothing."""

        pass

    def qcstep(self):
        """Dummy centroid position propagator which does nothing."""

        pass

    def step(self, step=None):
        """Dummy simulation time step which does nothing."""
        pass

    def pconstraints(self):
        pass


class NVEIntegrator(DummyIntegrator):
    """ Integrator object for constant energy simulations.

    Has the relevant conserved quantity and normal mode propagator for the
    constant energy ensemble. Note that a temperature of some kind must be
    defined so that the spring potential can be calculated.
    """
 
    def pstep(self, level=-1, alpha=1.0):
        """Velocity Verlet monemtum propagator."""
        
        if self.splitting == "aboba": dt = self.dt
        elif self.splitting == "obabo" or self.splitting == "baoab": dt = self.halfdt
        
        if level < 0: 
            self.beads.p += self.forces.f*(dt/alpha)
        else: 
            self.beads.p += self.forces.forces_mts(level)*(dt/alpha)
       
    def qcstep(self, alpha=1.0):
        """Velocity Verlet centroid position propagator."""
        if self.splitting == "obabo" or self.splitting == "baoab": dt = self.dt
        elif self.splitting == "aboba": dt = self.halfdt
        
        self.nm.qnm[0,:] += depstrip(self.nm.pnm)[0,:]/depstrip(self.beads.m3)[0]*dt/alpha
       
    def mtsprop(self, index, alpha):
        """ Recursive MTS step """
        nmtslevels = len(self.nmts)
        mk = self.nmts[index]  # mtslevels starts at level zero, where nmts should be 1 in most cases
        alpha *= mk
        
        if self.splitting == "obabo" or self.splitting == "baoab":  # fully equivalent
            for i in range(mk):  
                # propagate p for dt/2alpha with force at level index      
                self.pstep(index, alpha)
                self.pconstraints()
                
                if index == nmtslevels-1:
                # call Q propagation for dt/alpha at the inner step
                    self.qcstep(alpha)
                    self.nm.free_qstep() # this has been hard-wired to use the appropriate time step with depend magic
                else:
                    self.mtsprop(index+1, alpha)
     
                # propagate p for dt/2alpha
                self.pstep(index, alpha)
                self.pconstraints()
        elif self.splitting == "aboba":
            raise ValueError("Cannot use MTS with ABOBA splitting")
                            
                                
    def step(self, step=None):
        """Does one simulation time step."""
 
        if self.splitting == "obabo" or self.splitting == "baoab":
            # bias is applied at the outer loop too
            self.beads.p += depstrip(self.bias.f)*(self.halfdt)
            self.pconstraints()
        
            self.mtsprop(0,1.0)
 
            self.beads.p += depstrip(self.bias.f)*(self.halfdt)
            self.pconstraints()
        elif self.splitting == "aboba":   # cannot see how to MTS this in a meaningful way
            self.qcstep()
            self.nm.free_qstep() 
                    
            self.pstep()
            self.beads.p += depstrip(self.bias.f)*(self.dt)
            self.pconstraints()
                
            self.qcstep()
            self.nm.free_qstep()              

class NVTIntegrator(NVEIntegrator):
    """Integrator object for constant temperature simulations.

    Has the relevant conserved quantity and normal mode propagator for the
    constant temperature ensemble. Contains a thermostat object containing the
    algorithms to keep the temperature constant.

    Attributes:
        thermostat: A thermostat object to keep the temperature constant.
    """
 
    def pstep(self, level=0, alpha=1.0):
        """Velocity Verlet monemtum propagator."""
        
        # since this  is thermostatted, should use half dt for every splitting
        self.beads.p += self.forces.forces_mts(level)*(self.halfdt/alpha)
        if level == 0:  # bias is applied at the outer loop too
            self.beads.p += depstrip(self.bias.f)*(self.halfdt/alpha)

       
    def qcstep(self, alpha=1.0):
        """Velocity Verlet centroid position propagator."""
        if self.splitting == "obabo": dt = self.dt
        elif self.splitting == "aboba" or self.splitting == "baoab": dt = self.halfdt
        
        self.nm.qnm[0,:] += depstrip(self.nm.pnm)[0,:]/depstrip(self.beads.m3)[0]*dt/alpha
       
    def tstep(self):
        """Velocity Verlet thermostat step"""
        # the length of the thermostat step is controlled via depend objects
        self.thermostat.step()
        
    def mtsprop(self, index, alpha):
        """ Recursive MTS step """
        nmtslevels = len(self.nmts)
        mk = self.nmts[index]  # mtslevels starts at level zero, where nmts should be 1 in most cases
        alpha *= mk
        
        if self.splitting == "obabo":
            for i in range(mk):  
                # propagate p for dt/2alpha with force at level index      
                self.pstep(index, alpha)
                self.pconstraints()
                
                if index == nmtslevels-1:
                # call Q propagation for dt/alpha at the inner step
                    self.qcstep(alpha)
                    self.nm.free_qstep() # this has been hard-wired to use the appropriate time step with depend magic
                else:
                    self.mtsprop(index+1, alpha)
     
                # propagate p for dt/2alpha
                self.pstep(index, alpha)
                self.pconstraints()
        elif self.splitting == "baoab":
            for i in range(mk):  
                # propagate p for dt/2alpha with force at level index      
                self.pstep(index, alpha)
                self.pconstraints()
                
                if index == nmtslevels-1:
                    # call Q propagation for dt/2alpha at the inner step
                    self.qcstep(alpha)
                    self.nm.free_qstep()  # this has been hard-wired to use the appropriate time step with depend magic
                    self.tstep()
                    self.pconstraints()
                    self.qcstep(alpha)
                    self.nm.free_qstep()  # this has been hard-wired to use the appropriate time step with depend magic             
                else:
                    self.mtsprop(index+1, alpha)
     
                # propagate p for dt/2alpha
                self.pstep(index, alpha)
                self.pconstraints()
                
    def step(self, step=None):
        """Does one simulation time step."""
 
        if self.splitting == "obabo":
            # thermostat is applied at the outer loop
            self.tstep()
            self.pconstraints()
      
            self.mtsprop(0,1.0)
 
            self.tstep()
            self.pconstraints()
        elif self.splitting == "baoab":
            # bias is applied at the outer loop too            
            self.mtsprop(0,1.0)
            
        elif self.splitting == "aboba":   # cannot see how to MTS this in a meaningful way, so we do that manually
            self.qcstep()
            self.nm.free_qstep() 
                    
            self.pstep()
            self.pconstraints()

            self.tstep()
            self.pconstraints()

            self.pstep()
            self.pconstraints()
                
            self.qcstep()
            self.nm.free_qstep() 
            

class NPTIntegrator(NVTIntegrator):
    """Integrator object for constant pressure simulations.

    Has the relevant conserved quantity and normal mode propagator for the
    constant pressure ensemble. Contains a thermostat object containing the
    algorithms to keep the temperature constant, and a barostat to keep the
    pressure constant.
    """

    # should be enough to redefine these functions, and the step() from NVTIntegrator should do the trick
    def pstep(self, level=0, alpha=1.0):
        """Velocity Verlet monemtum propagator."""
        
        # since this  is thermostatted, should use half dt for every splitting
        self.barostat.pstep(level, 1.0/alpha)
        
       
    def qcstep(self, alpha=1.0):
        """Velocity Verlet centroid position propagator."""
        if self.splitting == "obabo": dt = 2.0
        elif self.splitting == "aboba" or self.splitting == "baoab": dt = 1.0
        
        self.barostat.qcstep(dtscale=dt/alpha)            
       
    def tstep(self):
        """Velocity Verlet thermostat step"""
        # the length of the thermostat step is controlled via depend objects
        self.thermostat.step()
        self.barostat.thermostat.step()
        

class NSTIntegrator(NVTIntegrator):
    """Ensemble object for constant pressure simulations.

    Has the relevant conserved quantity and normal mode propagator for the
    constant pressure ensemble. Contains a thermostat object containing the
    algorithms to keep the temperature constant, and a barostat to keep the
    pressure constant.

    Attributes:
    barostat: A barostat object to keep the pressure constant.

    Depend objects:
    econs: Conserved energy quantity. Depends on the bead and cell kinetic
    and potential energy, the spring potential energy, the heat
    transferred to the beads and cell thermostat, the temperature and
    the cell volume.
    pext: External pressure.
    """

    def step(self, step=None):
        """NST time step.

        Note that the barostat only propagates the centroid coordinates. If this
        approximation is made a centroid virial pressure and stress estimator can
        be defined, so this gives the best statistical convergence. This is
        allowed as the normal mode propagation is approximately unaffected
        by volume fluctuations as long as the system box is much larger than
        the radius of gyration of the ring polymers.
        """

        self.ttime = -time.time()
        if self.splitting == "obabo":
            self.thermostat.step()
            self.barostat.thermostat.step()
            self.pconstraints()
            self.barostat.pstep()
            self.pconstraints()

            self.barostat.qcstep(dtscale=2.0)
            self.nm.free_qstep()
            
            self.barostat.pstep()
            self.pconstraints()        
            self.barostat.thermostat.step()
            self.thermostat.step()
            self.pconstraints()
        elif self.splitting == "aboba":
            self.barostat.qcstep()
            self.nm.free_qstep()
            self.barostat.pstep()
            self.pconstraints()

            self.barostat.thermostat.step()
            self.thermostat.step()
            self.barostat.thermostat.step()
            self.pconstraints()
            
            self.barostat.pstep()
            self.pconstraints()
            
            self.barostat.qcstep()
            self.nm.free_qstep()
        elif self.splitting == "baoab":
            self.barostat.pstep()
            self.pconstraints()
            
            self.barostat.qcstep()
            self.nm.free_qstep()

            self.thermostat.step()
            self.pconstraints()
            
            self.barostat.qcstep()
            self.nm.free_qstep()
        
            self.barostat.pstep()
            self.pconstraints()
            
        self.ttime += time.time()

class SCIntegrator(NVTIntegrator):
    """Integrator object for constant temperature simulations.

    Has the relevant conserved quantity and normal mode propagator for the
    constant temperature ensemble. Contains a thermostat object containing the
    algorithms to keep the temperature constant.

    Attributes:
        thermostat: A thermostat object to keep the temperature constant.

    Depend objects:  
        econs: Conserved energy quantity. Depends on the bead kinetic and
            potential energy, the spring potential energy and the heat
            transferred to the thermostat.
    """

    def bind(self, mover):
        """Binds ensemble beads, cell, bforce, bbias and prng to the dynamics.

        This takes a beads object, a cell object, a forcefield object and a
        random number generator object and makes them members of the ensemble.
        It also then creates the objects that will hold the data needed in the
        ensemble algorithms and the dependency network. Note that the conserved
        quantity is defined in the init, but as each ensemble has a different
        conserved quantity the dependencies are defined in bind.

        Args:
        beads: The beads object from whcih the bead positions are taken.
        nm: A normal modes object used to do the normal modes transformation.
        cell: The cell object from which the system box is taken.
        bforce: The forcefield object from which the force and virial are
            taken.
        prng: The random number generator object which controls random number
            generation.
        """
      
        super(SCIntegrator,self).bind(mover)
        self.ensemble.add_econs(dget(self.forces, "potsc"))

    def pstep(self):                                                                     
        """Velocity Verlet momenta propagator."""
                              
                                                
        # also include the baseline Tr2SC correction (the 2/3 & 4/3 V bit)
        self.beads.p += depstrip(self.forces.f + self.coeffsc*self.forces.f)*self.halfdt/self.nmts
        # also adds the bias force (TODO!!!)
        # self.beads.p += depstrip(self.bias.f)*(self.dt*0.5)
                                                                                        
    def pscstep(self):                                                                     
        """Velocity Verlet Suzuki-Chin momenta propagator."""

        # also adds the force assiciated with SuzukiChin correction (only the |f^2| term, so we remove the Tr2SC correction)
        self.beads.p += depstrip(self.forces.fsc - self.coeffsc*self.forces.f)*self.halfdt

    def qcstep(self):
        """Velocity Verlet centroid position propagator."""
                                 
        if self.splitting == "aboba" or self.splitting == "baoab": dt = self.halfdt
        elif self.splitting == "obabo": dt = self.dt      
                                                                                         
        self.nm.qnm[0,:] += depstrip(self.nm.pnm)[0,:]/depstrip(self.beads.m3)[0]*dt/self.nmts

    def step(self, step=None):
        """Does one simulation time step."""


        if self.splitting == "obabo":
            self.thermostat.step()
            self.pconstraints()
            
            self.pscstep()

            for i in range(self.nmts):
                self.pstep()
                self.pconstraints()
          
                self.qcstep()
                self.nm.free_qstep()
           
                self.pstep()
          
            self.pscstep()
            self.pconstraints()

            self.thermostat.step()
            self.pconstraints()
        elif self.splitting == "baoab":
            self.pscstep()
            for i in range(self.nmts):
                self.pstep()
                self.pconstraints()
          
                self.qcstep()
                self.nm.free_qstep()
                self.thermostat.step()
                self.pconstraints()
                self.qcstep()
                self.nm.free_qstep()
                
                self.pstep()
          
            self.pscstep()
            self.pconstraints()
