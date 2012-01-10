"""
Contains the classes that connect the driver to the python code.

Communicates with the driver code, obtaining the force, virial and potential.
Deals with creating the jobs that will be sent to the driver, and 
returning the results to the python code.

Classes:
   RestartForce: Deals with creating the ForceField object from a file, and
      writing the checkpoints.
   ForceField: Base forcefield class with the generic methods and attributes.
   FFSocket: Deals with a single replica of the system
   ForceBeads: Deals with the parallelization of the force calculation for
      a PI simulation.

Exceptions:

Functions:

"""

__all__ = ['RestartForce', 'ForceField', 'ForceBeads', 'FFSocket']

import numpy as np
import math, time, threading
from utils.depend import *
from utils.io import io_system
from driver.interface import Interface, RestartInterface
from utils.restart import *

class RestartForce(Restart):
   """Forcefield restart class.

      Handles generating the appropriate forcefield class from the xml
      input file, and generating the xml checkpoint tags and data from an 
      instance of the object.

      Attributes:
         attribs: A dictionary giving the attributes. Of the form 
            {"name": (type, (data type, default value))}
         fields: A dictionary giving the fields. Of the form 
            {"name": (type, (data type, default value))}

      Attribs:
         type: A string indicating the type being used. 'socket' is currently
            the only available option.

      Fields:
         interface: A restart interface instance.
         parameters: A dictionary of the parameters used by the driver. Of the
            form {"name": value}.
      """

   attribs = { "type" : (RestartValue,(str,"socket")) }
   fields =  { "interface" : (RestartInterface,()), "parameters" : (RestartValue, (dict,None)) }
   
   def store(self, force):
      """Takes a ForceField instance and stores a minimal representation of it.

      Args:
         force: A forcefield object.
      """

      if (type(force) is FFSocket):  
         self.type.store("socket")
         self.interface.store(force.socket)
         self.parameters.store(force.pars)
      else: self.type.store("unknown")
         

   def fetch(self):
      """Creates a forcefield object.

      Returns:
         A forcefield object of the appropriate type and with the appropriate
         interface given the attributes of the RestartForce object.
      """

      if self.type.fetch().upper() == "SOCKET": 
         force=FFSocket(pars=self.parameters.fetch(), interface=self.interface.fetch())
      else : force=ForceField()
      return force
      

class ForceField(dobject):
   """Base forcefield class.

   Gives the standard methods and quantities needed in all the forcefield
   classes.

   Attributes: 
      ufv: A list of the form [pot, f, vir].
      pot = The potential energy of the system.
      f: An array containing all the components of the force.
      fx: A slice of f containing only the x components of the forces.
      fy: A slice of f containing only the y components of the forces.
      fz: A slice of f containing only the z components of the forces.
      vir: An array containing the components of the virial tensor, 
         not divided by the volume.
   """

   def __init__(self):
      """Initialises ForceField."""

      dset(self,"ufv", depend_value(name="ufv", func=self.get_all))
      
   def copy(self):
      """Creates a deep copy without the bound objects.

      Returns:
         A ForceField object without atoms or cell attributes.
      """

      return type(self)()
      
   def bind(self, atoms, cell):
      """Binds atoms and cell to the forcefield.

      This takes an atoms object and a cell object and makes them members of
      the forcefield. It also then creates the objects that will hold the data
      that the driver returns and the dependency network.

      Args:
         atoms: The atoms from which the positions are taken.
         cell: The cell from which the cell box is taken.
      """

      self.atoms = atoms
      self.cell = cell
      dget(self,"ufv").add_dependency(dget(self.atoms,"q"))
      dget(self,"ufv").add_dependency(dget(self.cell,"h"))      
      dset(self,"pot",depend_value(name="pot", func=self.get_pot, dependencies=[dget(self,"ufv")] )  )
      
      fbase=np.zeros(atoms.natoms*3, float)
      dset(self,"f", depend_array(name="f", value=fbase, func=self.get_f, dependencies=[dget(self,"ufv")]) )
      # a bit messy, but we don't want to trigger quite yet the force calculation routine 
      dset(self,"fx", depend_array(name="fx", value=fbase[0:3*atoms.natoms:3]));
      dset(self,"fy", depend_array(name="fy", value=fbase[1:3*atoms.natoms:3]));
      dset(self,"fz", depend_array(name="fz", value=fbase[2:3*atoms.natoms:3]));
      depcopy(self,"f", self,"fx");      depcopy(self,"f", self,"fy");      depcopy(self,"f", self,"fz");
      dset(self,"vir", depend_array(name="vir", value=np.zeros((3,3),float),func=self.get_vir, dependencies=[dget(self,"ufv")] ) )

   def queue(self):
      """Dummy queueing method."""
      
      pass
   
   def get_all(self):
      """Dummy driver routine.

      Returns:
         A list where the potential, force and virial have had all 
         components set to zero.
      """

      return [0.0, numpy.zeros(3*self.atoms.natoms), numpy.zeros((3,3),float)]

   def get_pot(self):
      """Calls get_all routine of forcefield to update potential."""

      [pot, f, vir] = self.ufv
      return pot

   def get_f(self):
      """Calls get_all routine of forcefield to update force."""

      [pot, f, vir] = self.ufv
      return f

   def get_vir(self):
      """Calls get_all routine of forcefield to update virial."""

      [pot, f, vir] = self.ufv
      vir[1,0]=0.0; vir[2,0:2]=0.0;
      return vir


class ForceBeads(dobject):
   """PI forcefield object.

   Collects many forcefield instances and parallelizes getting the forces 
   in a PIMD environment. Deals with splitting the bead representation into
   separate replicas, and collecting the data from each replica.

   Attributes:
      natoms: Number of atoms.
      nbeads: Number of beads.
      _forces: A list containing all the force objects for each system replica.
      f: An array containing the components of the force.
      pots: A list containing the potential energy for each system replica.
      virs: A list containing the virial tensor for each system replica.
      pot: The appropriate estimator for the potential energy.
      vir: The appropriate estimator for the virial tensor.
      fnm: An array containing the components of the force in the normal mode
         representation.
   """

   def __init__(self, beads=None, cell=None, force=None):
      """Initialises ForceBeads.

      Args:
         beads: Optional beads object, to be bound to the forcefield.
         cell: Optional cell object, to be bound to the forcefield.
         force: Optional force object, to be bound to the forcefield. This
            should be a FFSocket or equivalent, so that it can be copied for
            each replica of the system.
      """

      if not (beads is None or cell is None or force is None): 
         self.bind(beads, cell, force)
      else: 
         pass   

   def bind(self, beads, cell, force):
      """Binds beads, cell and force to the forcefield.

      Takes the beads, cell objects and makes them members of the forcefield.
      Also takes the force object and copies it once for each replica of the
      system, then binds each replica to one of the copies so that the force
      calculation can be parallelized. Creates the objects that will 
      hold the data that the driver returns and the dependency network. 

      Args:
         beads: Beads object, to be bound to the forcefield.
         cell: Cell object, to be bound to the forcefield.
         force: Force object, to be bound to the forcefield. This
            should be a FFSocket or equivalent, so that it can be copied for
            each replica of the system.
      """

      self.natoms=beads.natoms
      self.nbeads=beads.nbeads

      self._forces=[];
      for b in range(self.nbeads):
         newf=force.copy()
         newf.bind(beads[b], cell)
         newf.blocking=False
         self._forces.append(newf)
               
      dset(self,"f",depend_array(name="f",value=np.zeros((self.nbeads,3*self.natoms), float), func=self.f_gather,     
          dependencies=[dget(self._forces[b],"f")  for b in range(self.nbeads)] ) )
      dset(self,"pots",depend_array(name="pots", value=np.zeros(self.nbeads,float), func=self.pot_gather,     
          dependencies=[dget(self._forces[b],"pot")  for b in range(self.nbeads)] ) )
      dset(self,"virs",depend_array(name="virs", value=np.zeros((self.nbeads,3,3),float), func=self.vir_gather,     
          dependencies=[dget(self._forces[b],"vir")  for b in range(self.nbeads)] ) )          
      dset(self,"pot",depend_value(name="pot", func=self.pot,     
          dependencies=[dget(self,"pots")] ) )
      dset(self,"vir",depend_value(name="vir", func=self.vir,
          dependencies=[dget(self,"virs")] ) )

      dset(self,"fnm",depend_array(name="fnm",value=np.zeros((self.nbeads,3*self.natoms), float), func=self.b2nm_f, dependencies=[dget(self,"f")] ) )
      self.Cb2nm=beads.Cb2nm
      
   def b2nm_f(self): return np.dot(self.Cb2nm,depstrip(self.f))

   def _getbead(self, b, newf):
      newf[b]=self._forces[b].f
      return

   def queue(self):   
      for b in range(self.nbeads): self._forces[b].queue()

   def pot_gather(self): 
      self.queue()
      return np.array([b.pot for b in self._forces], float)

   def vir_gather(self): 
      self.queue()
      return np.array([b.vir for b in self._forces], float)
      
   def pot(self): return self.pots.sum()

   def vir(self): return self.virs.sum()

   def f_gather(self): 
      start=time.time()
      newf=np.zeros((self.nbeads,3*self.natoms),float)
      
      self.queue()

      #serial
      for b in range(self.nbeads): newf[b]=self._forces[b].f
      # threaded      

      return newf


LATENCY=5e-3
class FFSocket(ForceField):

   def __init__(self, pars={}, interface=None):
      super(FFSocket,self).__init__() 
      if interface is None:
         self.socket=Interface()
      else:
         self.socket=interface
      self.pars=pars     
      self.request=None
      
      self.timer=0.0
      self.twall=0.0
      self.ncall=0
      
   def copy(self):    # creates a deep copy with everything but the bound bits 
      return type(self)(self.pars, self.socket)

   def get_all(self):
      self.timer-= time.clock()
      self.twall-=time.time()
      
      if self.request is None: self.request=self.socket.queue(self.atoms, self.cell, self.pars)
      while self.request["status"] != "Done": time.sleep(self.socket.latency)
      self.socket.release(self.request)
      result=self.request["result"]
      self.request=None
      
      self.ncall+=1
      self.timer+=time.clock()
      self.twall+=time.time()      
      return result
      
   def queue(self):   
      if self.request is None and dget(self,"ufv").tainted():
         self.request=self.socket.queue(self.atoms, self.cell, self.pars)
