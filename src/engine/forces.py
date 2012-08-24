"""Contains the classes that connect the driver to the python code.

Communicates with the driver code, obtaining the force, virial and potential.
Deals with creating the jobs that will be sent to the driver, and
returning the results to the python code.

Classes:
   ForceField: Base forcefield class with the generic methods and attributes.
   FFSocket: Deals with a single replica of the system
   ForceBeads: Deals with the parallelization of the force calculation for
      a PI simulation.
"""

__all__ = ['ForceField', 'ForceBeads', 'Forces', 'FFSocket']

import numpy as np
import math, time
from utils.depend import *
from utils.nmtransform import nm_rescale
from driver.interface import Interface
from beads import Beads


class ForceField(dobject):
   """Base forcefield class.

   Gives the standard methods and quantities needed in all the forcefield
   classes.

   Attributes:
      atoms: An Atoms object containing all the atom positions.
      cell: A Cell object containing the system box.
      softexit: A function to help make sure the printed restart file is
         consistent.

   Depend objects:
      ufv: A list of the form [pot, f, vir]. These quantities are calculated
         all at one time by the driver, so are collected together. Each separate
         object is then taken from the list. Depends on the atom positions and
         the system box.
      pot: A float giving the potential energy of the system. Depends on ufv.
      f: An array containing all the components of the force. Depends on ufv.
      fx: A slice of f containing only the x components of the forces.
      fy: A slice of f containing only the y components of the forces.
      fz: A slice of f containing only the z components of the forces.
      vir: An array containing the components of the virial tensor in upper
         triangular form, not divided by the volume. Depends on ufv.
   """

   def __init__(self):
      """Initialises ForceField."""

      # ufv is a list [ u, f, vir ]  which stores the results of the force calculation
      dset(self,"ufv", depend_value(name="ufv", func=self.get_all))

   def copy(self):
      """Creates a deep copy without the bound objects.

      Used in ForceBeads to create a ForceField for each replica of the system.

      Returns:
         A ForceField object without atoms or cell attributes.
      """

      return type(self)(self.nbeads, self.weight)

   def bind(self, atoms, cell, softexit=None):
      """Binds atoms and cell to the forcefield.

      This takes an atoms object and a cell object and makes them members of
      the forcefield. It also then creates the objects that will hold the data
      that the driver returns and the dependency network.

      Args:
         atoms: The Atoms object from which the atom positions are taken.
         cell: The Cell object from which the system box is taken.
         softexit: A function to help make sure the printed restart file is
            consistent.
      """

      # stores a reference to the atoms and cell we are computing forces for
      self.atoms = atoms
      self.cell = cell
      self.softexit = softexit

      # ufv depends on the atomic positions and on the cell
      dget(self,"ufv").add_dependency(dget(self.atoms,"q"))
      dget(self,"ufv").add_dependency(dget(self.cell,"h"))

      # potential and virial are to be extracted very simply from ufv
      dset(self,"pot",
         depend_value(name="pot", func=self.get_pot,
            dependencies=[dget(self,"ufv")]))

      dset(self,"vir",
         depend_array(name="vir", value=np.zeros((3,3),float),func=self.get_vir,
            dependencies=[dget(self,"ufv")]))


       # NB: the force requires a bit more work, to define shortcuts to xyz slices
       # without calculating the force at this point.
      fbase = np.zeros(atoms.natoms*3, float)
      dset(self,"f",
         depend_array(name="f", value=fbase, func=self.get_f,
             dependencies=[dget(self,"ufv")]))

      dset(self,"fx", depend_array(name="fx", value=fbase[0:3*atoms.natoms:3]))
      dset(self,"fy", depend_array(name="fy", value=fbase[1:3*atoms.natoms:3]))
      dset(self,"fz", depend_array(name="fz", value=fbase[2:3*atoms.natoms:3]))
      depcopy(self,"f", self,"fx")
      depcopy(self,"f", self,"fy")
      depcopy(self,"f", self,"fz")



   def queue(self):
      """Dummy queueing method."""

      pass

   def stop(self):
      """Dummy queueing method."""

      pass

   def run(self):
      """Dummy queueing method."""

      pass

   def get_all(self):
      """Dummy driver routine.

      Returns:
         A list of the form [potential, force, virial] where the potential
         and all components of the force and virial have been set to zero.
      """

      return [0.0, np.zeros(3*self.atoms.natoms), np.zeros((3,3),float)]

   def get_pot(self):
      """Calls get_all routine of forcefield to update potential.

      Returns:
         Potential energy.
      """

      return self.ufv[0]

   def get_f(self):
      """Calls get_all routine of forcefield to update force.

      Returns:
         An array containing all the components of the force.
      """

      return self.ufv[1]

   def get_vir(self):
      """Calls get_all routine of forcefield to update virial.

      Returns:
         An array containing the virial in upper triangular form, not divided
         by the volume.
      """

      vir = self.ufv[2]
      vir[1,0] = 0.0
      vir[2,0:2] = 0.0
      return vir


class FFSocket(ForceField):
   """Interface between the PIMD code and the socket for a single replica.

   Deals with an individual replica of the system, obtaining the potential
   force and virial appropriate to this system. Deals with the distribution of
   jobs to the interface.

   Attributes:
      parameters: A dictionary of the parameters used by the driver. Of the
         form {'name': value}.
      socket: The interface object which contains the socket through which
         communication between the forcefield and the driver is done.
      request: During the force calculation step this holds a dictionary
         containing the relevant data for determining the progress of the step.
         Of the form {'atoms': atoms, 'cell': cell, 'pars': parameters,
                      'status': status, 'result': result, 'id': bead id,
                      'start': starting time}.
   """

   def __init__(self, pars=None, interface=None):
      """Initialises FFSocket.

      Args:
         pars: Optional dictionary, giving the parameters needed by the driver.
         interface: Optional Interface object, which contains the socket.
      """

      # a socket to the communication library is created or linked
      super(FFSocket,self).__init__()
      if interface is None:
         self.socket = Interface()
      else:
         self.socket = interface

      if pars is None:
         self.pars = {}
      else:
         self.pars = pars
      self.request = None

   def bind(self, atoms, cell, softexit=None):
      """Pass on the binding request from ForceBeads.

      Also makes sure to set the socket's softexit.

      Args:
         atoms: Atoms object from which the bead positions are taken.
         cell: Cell object from which the system box is taken.
         softexit: A function to help make sure the printed restart file is
            consistent.
      """

      super(FFSocket,self).bind(atoms, cell, softexit)
      if not self.softexit is None:
         self.socket.softexit = self.softexit

   def copy(self):
      """Creates a deep copy without the bound objects.

      Used in ForceBeads to create a FFSocket for each replica of the system.

      Returns:
         A FFSocket object without atoms or cell attributes.
      """

      # does not copy the bound objects (i.e., the returned forcefield must be bound before use)
      return type(self)(self.pars, self.socket)

   def get_all(self):
      """Driver routine.

      When one of the force, potential or virial are called, this sends the
      atoms and cell to the driver through the interface, requesting that the
      driver does the calculation. This then waits until the driver is finished,
      and then returns the ufv list.

      Returns:
         A list of the form [potential, force, virial].
      """

      # this is converting the distribution library requests into [ u, f, v ]  lists
      if self.request is None:
         self.request = self.socket.queue(self.atoms, self.cell, pars=self.pars, reqid=-1)
      while self.request["status"] != "Done":
         if self.request["status"] == "Exit":
            break
         time.sleep(self.socket.latency)
      if self.request["status"] == "Exit":
         print " @Force:   Soft exit request."
         if not self.softexit is None:
            self.softexit()

      # data has been collected, so the request can be released and a slot freed up for new calculations
      self.socket.release(self.request)
      result = self.request["result"]
      self.request = None

      return result

   def queue(self, reqid=-1):
      """Sends the job to the interface queue directly.

      Allows the ForceBeads object to ask for the ufv list of each replica
      directly without going through the get_all function. This allows
      all the jobs to be sent at once, allowing them to be parallelized.

      Args:
         reqid: An optional integer that indentifies requests of the same type,
            e.g. the bead index.
      """

      if self.request is None and dget(self,"ufv").tainted():
         self.request = self.socket.queue(self.atoms, self.cell, pars=self.pars, reqid=reqid)

   def run(self):
      if not self.socket.started(): self.socket.start_thread()

   def stop(self):
      if self.socket.started(): self.socket.end_thread()


class ForceBeads(dobject):
   """Class that gathers all the forces together.

   Collects many forcefield instances and parallelizes getting the forces
   in a PIMD environment. Deals with splitting the bead representation into
   separate replicas, and collecting the data from each replica.

   Attributes:
      natoms: An integer giving the number of atoms.
      nbeads: An integer giving the number of beads.
      _forces: A list containing all the force objects for each system replica.
      softexit: A function to help make sure the printed restart file is
         consistent.

   Depend objects:
      f: An array containing the components of the force. Depends on each
         replica's ufv list.
      pots: A list containing the potential energy for each system replica.
         Depends on each replica's ufv list.
      virs: A list containing the virial tensor for each system replica.
         Depends on each replica's ufv list.
      pot: The sum of the potential energy of the replicas.
      vir: The sum of the virial tensor of the replicas.
      fnm: An array containing the components of the force in the normal mode
         representation.
   """

   def __init__(self, model, nbeads=0, weight=1.0):
      self.f_model = model
      self.nbeads = nbeads
      self.weight = weight

   def copy(self):
      """Creates a deep copy without the bound objects.

      Used in ForceBeads to create a FFSocket for each replica of the system.

      Returns:
         A FFSocket object without atoms or cell attributes.
      """

      # does not copy the bound objects (i.e., the returned forcefield must be bound before use)
      return type(self)(self.f_model, self.nbeads, self.weight)


   def bind(self, beads, cell, softexit=None):
      """Binds beads, cell and force to the forcefield.

      Takes the beads, cell objects and makes them members of the forcefield.
      Also takes the force object and copies it once for each replica of the
      system, then binds each replica to one of the copies so that the force
      calculation can be parallelized. Creates the objects that will
      hold the data that the driver returns and the dependency network.

      Args:
         beads: Beads object from which the bead positions are taken.
         cell: Cell object from which the system box is taken.
         force: Force object, to be bound to the forcefield. This
            should be a FFSocket or equivalent, so that it can be copied for
            each replica of the system.
         softexit: A function to help make sure the printed restart file is
            consistent.
      """

      # stores a copy of the number of atoms and of beads !TODO! make them read-only properties
      self.natoms = beads.natoms
      self.softexit = softexit
      if (self.nbeads != beads.nbeads) : raise ValueError("Binding together a Beads and a ForceBeads objects with different numbers of beads")

      # creates an array of force objects, which are bound to the beads and the cell
      self._forces = [];
      for b in range(self.nbeads):
         new_force = self.f_model.copy()
         new_force.bind(beads[b], cell, softexit=self.softexit)
         self._forces.append(new_force)

      # f is a big array which assembles the forces on individual beads
      dset(self,"f",
         depend_array(name="f",value=np.zeros((self.nbeads,3*self.natoms)),
            func=self.f_gather,
               dependencies=[dget(self._forces[b],"f") for b in range(self.nbeads)]))

      # collection of pots and virs from individual beads
      dset(self,"pots",
         depend_array(name="pots", value=np.zeros(self.nbeads,float),
            func=self.pot_gather,
               dependencies=[dget(self._forces[b],"pot") for b in range(self.nbeads)]))
      dset(self,"virs",
         depend_array(name="virs", value=np.zeros((self.nbeads,3,3),float),
            func=self.vir_gather,
               dependencies=[dget(self._forces[b],"vir") for b in range(self.nbeads)]))

      # total potential and total virial
      dset(self,"pot",
         depend_value(name="pot", func=self.get_pot,
            dependencies=[dget(self,"pots")]))
      dset(self,"vir",
         depend_value(name="vir", func=self.get_vir,
            dependencies=[dget(self,"virs")]))

   def run(self):
      for b in range(self.nbeads):
         self._forces[b].run()

   def stop(self):
      for b in range(self.nbeads):
         self._forces[b].stop()

   def queue(self):
      """Submits all the required force calculations to the interface."""

      # this should be called in functions which access u,v,f for ALL the beads,
      # before accessing them. it is basically pre-queueing so that the
      # distributed-computing magic can work
      for b in range(self.nbeads):
         self._forces[b].queue(reqid=b)


   def pot_gather(self):
      """Obtains the potential energy for each replica.

      Returns:
         A list of the potential energy of each replica of the system.
      """

      self.queue()
      return np.array([b.pot for b in self._forces], float)

   def vir_gather(self):
      """Obtains the virial for each replica.

      Returns:
         A list of the virial of each replica of the system.
      """

      self.queue()
      return np.array([b.vir for b in self._forces], float)


   def f_gather(self):
      """Obtains the global force vector.

      Returns:
         An array with all the components of the force. Row i gives the force
         array for replica i of the system.
      """

      newf = np.zeros((self.nbeads,3*self.natoms),float)

      self.queue()
      for b in range(self.nbeads):
         newf[b] = self._forces[b].f

      #serial
#      for b in range(self.nbeads): newf[b]=self._forces[b].f
      # threaded
#      bthreads=[]
#      print "starting threads"
#      for b in range(self.nbeads):
#         thread=threading.Thread(target=self._getbead, args=(b,newf,))
#         thread.start()
#         bthreads.append(thread)

#      print "waiting threads"
#      for b in range(self.nbeads): bthreads[b].join()
#      print "threads joined in"

      return newf

   def get_pot(self):
      """Sums the potential of each replica.

      Used to check energy conservation. Not the actual system potential.

      Returns:
         Potential energy sum.
      """

      return self.pots.sum()

   def get_vir(self):
      """Sums the virial of each replica.

      Not the actual system virial.

      Returns:
          Virial sum.
      """

      vir = np.zeros((3,3))
      for v in self.virs:
         vir += v
      return vir

   def __len__(self):
      """Length function.

      This is called whenever the standard function len(forcebeads) is used.

      Returns:
         The number of beads.
      """

      return self.nbeads

   def __getitem__(self,index):
      """Overwrites standard getting function.

      This is called whenever the standard function forcebeads[index] is used.
      Returns the force on bead index.

      Args:
         index: The index of the replica of the system to be accessed.

      Returns:
         The replica of the system given by the index.
      """

      return self._forces[index]


class Forces(dobject):

   def bind(self, beads, cell, flist, softexit=None):

      self.natoms = beads.natoms
      self.nbeads = beads.nbeads
      self.nforces = len(flist)
      self.softexit = softexit


      # flist should be a list of tuples containing ( "name", forcebeads)
      self.mforces=[]
      self.mbeads=[]
      self.mweights=[]
      self.mrpc=[]


      # a "function factory" to generate functions to automatically update contracted paths
      def make_rpc(rpc, beads):
         return lambda: rpc.b1tob2(depstrip(beads.q))


      # creates new force objects, possibly acting on contracted path representations
      for ( ftype, fbeads) in flist:

         # creates an automatically-updated contracted beads object
         newb=fbeads.nbeads
         newforce=fbeads.copy()
         newweight=fbeads.weight

         # if the number of beads for this force component is unspecified, assume full force evaluation
         if newb == 0: newb = beads.nbeads

         newbeads=Beads(beads.natoms, newb)
         newrpc=nm_rescale(beads.nbeads, newb)

         newf=make_rpc(newrpc, beads)
         dget(newbeads,"q")._func=newf
         for b in newbeads: dget(b,"q")._func=newf      # must update also indirect access to the beads coordinates
         dget(beads,"q").add_dependant(dget(newbeads,"q"))        # makes newbeads.q depend from beads.q

         #now we create a new forcebeads which is bound to newbeads!
         newforce.bind(newbeads, cell, softexit)

         self.mweights.append(newweight)
         self.mbeads.append(newbeads)
         self.mforces.append(newforce)
         self.mrpc.append(newrpc)

      #now must expose an interface that gives overall forces
      dset(self,"f",
         depend_array(name="f",value=np.zeros((self.nbeads,3*self.natoms)),
            func=self.f_combine,dependencies=[dget(ff, "f") for ff in self.mforces]
         ) )

      # collection of pots and virs from individual ff objects
      dset(self,"pots",
         depend_array(name="pots", value=np.zeros(self.nbeads,float),
            func=self.pot_combine,dependencies=[dget(ff, "pots") for ff in self.mforces]) )

      # must take care of the virials!
      dset(self,"virs",
         depend_array(name="virs", value=np.zeros((self.nbeads,3,3),float),
         func=self.vir_combine, dependencies=[dget(ff, "virs") for ff in self.mforces]) )

      # total potential and total virial
      dset(self,"pot",
         depend_value(name="pot", func=(lambda: self.pots.sum()),
            dependencies=[dget(self,"pots")]))
      dset(self,"vir",
         depend_value(name="vir", func=(lambda: self.virs.sum()),
            dependencies=[dget(self,"virs")]))

   def run(self):
      for ff in self.mforces: ff.run()

   def stop(self):
      for ff in self.mforces: ff.stop()

   def queue(self):
      for ff in self.mforces: ff.queue()

   def f_combine(self):

      self.queue()
      rf=np.zeros((self.nbeads,3*self.natoms),float)
      for k in range(self.nforces):
         rf+=self.mweights[k]*self.mrpc[k].b2tob1(self.mforces[k].f)   # "expand" to the total number of beads the forces from the contracted one
      return rf

   def pot_combine(self):

      self.queue()
      rp=np.zeros(self.nbeads,float)
      for k in range(self.nforces):
         rp+=self.mweights[k]*self.mrpc[k].b2tob1(self.mforces[k].pots)   # "expand" to the total number of beads the potentials from the contracted one
      return rp

   def vir_combine(self):

      self.queue()

      rp=np.zeros((self.nbeads,3,3),float)
      for k in range(self.nforces):
         virs=depstrip(self.mforces[k].virs)
         # "expand" to the total number of beads the virials from the contracted one, element by element
         for i in range(3):
            for j in range(3):
               rp[:,i,j]+=self.mweights[k]*self.mrpc[k].b2tob1(virs[:,i,j])
      return rp
