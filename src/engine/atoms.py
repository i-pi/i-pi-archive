import numpy as np
from utils.depend import *
from utils.restart import *
import utils.io.io_pdb
from utils import units

class Atom(dobject):
   """Represent an atom, with position, velocity, mass and related properties.

   This is actually only an interface to the Atoms class, i.e. only stores
   views of the large arrays which contain all the coordinates. 

   Attributes:
      p: The three components of the momentum of the atom.
      q: The three components of the position of the atom.
      m: The mass of the atom.
      name: The name of the atom.
      m3: An array of 3 elements with each element being the mass of the atom.
         Used when each degree of freedom needs to be divided by the mass.
      kin: The kinetic energy of the atom.
      kstress: The contribution of the atom to the kinetic stress tensor.
   """
            
   def __init__(self, system, index):
      """Initialises Atom.

      Args:
         system: An Atoms object containing the required atom.
         index: An integer giving the index of the required atom in the atoms
            list. Note that indices start from 0.
      """

      dset(self,"p",system.p[3*index:3*index+3])
      dset(self,"q",system.q[3*index:3*index+3])
      dset(self,"m",system.m[index:index+1])
      dset(self,"name",system.names[index:index+1])
      dset(self,"m3",system.m3[3*index:3*index+3])
      
   @property
   def kin(self):
      """Calculates the contribution of the atom to the kinetic energy."""

      return np.dot(self.p,self.p)/(2.0*self.m)

   @property
   def kstress(self):
      """Calculates the contribution of the atom to the kinetic stress 
      tensor.
      """

      p = depstrip(self.p)
      ks = numpy.zeros((3,3),float)
      for i in range(3):
         for j in range(i,3):
            ks[i,j] = p[i]*p[j]            
      return ks/self.m


class Atoms(dobject):
   """Storage for the atoms' positions, masses and velocities. 

   Everything is stored as 3*n sized contiguous arrays,
   and a convenience-access is provided through a list of Atom objects.

   Attributes:
      natoms: An integer giving the number of atoms.
      p: An array giving the components of the atom positions.
      q: An array giving the components of the atom momenta.
      m: An array giving the atom masses.
      names: An array giving the atom names.
      m3: An array of 3*n elements where each element of m has been copied 
         three times. Used when each degree of freedom needs to be divided 
         by the mass.
      M: The total mass of all the atoms.
      kin: The total kinetic energy of the atoms.
      kstress: The contribution of the atoms to the kinetic stress tensor.
      qx: An array giving the x components of the positions.
      qy: An array giving the y components of the positions.
      qz: An array giving the z components of the positions.
      px: An array giving the x components of the momenta.
      py: An array giving the y components of the momenta.
      pz: An array giving the z components of the momenta.
   """
            

   def __init__(self, natoms, _prebind=None):
      """Initialises Atoms.

      Each replica and the centroid coordinate are all held as Atoms objects, 
      and so slices of the global position and momentum arrays must be used in
      the initialisation so that they always agree with each other.

      Args:
         natoms: An integer giving the number of atoms.
         _prebind: An optional tuple of four elements; a depend_array of length
            3*natoms for the positions, another for the momenta, a depend_array
            of length natoms for the masses and another for the names.
      """

      self.natoms = natoms
      
      if _prebind is None:
         dset(self,"q",depend_array(name="q",value=np.zeros(3*natoms, float)) ) 
         dset(self,"p",depend_array(name="p",value=np.zeros(3*natoms, float)) )
         dset(self,"m",depend_array(name="m",value=np.zeros(natoms, float)) )
         dset(self,"names",depend_array(name="names",value=np.zeros(natoms, np.dtype('|S6'))) )         
      else:
         dset(self,"q",_prebind[0]) 
         dset(self,"p",_prebind[1]) 
         dset(self,"m",_prebind[2])
         dset(self,"names",_prebind[3])
 
      dset(self,"px",self.p[0:3*natoms:3],name="px")
      dset(self,"py",self.p[1:3*natoms:3],name="py")
      dset(self,"pz",self.p[2:3*natoms:3],name="pz")
      dset(self,"qx",self.q[0:3*natoms:3],name="qx")
      dset(self,"qy",self.q[1:3*natoms:3],name="qy")
      dset(self,"qz",self.q[2:3*natoms:3],name="qz")      
      
      dset(self,"m3",depend_array(name="m3",value=np.zeros(3*natoms, float),func=self.mtom3, dependencies=[dget(self,"m")]))

      dset(self,"M",depend_value(name="M",func=self.get_msum,dependencies=[dget(self,"m")]) )      
      dset(self,"kin",depend_value(name="kin",func=self.get_kin,dependencies=[dget(self,"p"),dget(self,"m3")]) )
      dset(self,"kstress",depend_value(name="kstress",func=self.get_kstress,dependencies=[dget(self,"p"),dget(self,"m")]) )
   
   def copy(self):
      """Creates a new Atoms object.

      Returns:
         An Atoms object with the same q, p, m and names arrays as the original.
      """

      newat = Atoms(self.natoms)
      newat.q[:] = self.q
      newat.p[:] = self.p
      newat.m[:] = self.m
      newat.names[:] = self.names
      return newat
      
   def __len__(self):
      """Length function.

      This is called whenever the standard function len(atoms) is used. Returns
      natoms.
      """

      return self.natoms

   def __getitem__(self,index):
      """Overwrites standard getting function.

      This is called whenever the standard function atoms[index] is used. 
      Returns an Atom object with the appropriate position and momenta arrays.
      Note that they are dynamically generated each time an Atom needs to be
      accessed, as this reduces the number of depend objects that need to be
      held at any one time.

      Args:
         index: The atom that needs to be accessed.
      """

      return Atom(self,index)

   def __setitem__(self,index,value):
      """Overwrites standard setting function.

      This is called whenever the standard function atoms[index]=value is used. 
      Changes the position and momenta of the appropriate slice of the global
      position and momentum arrays to those given by value.
      Note that they are dynamically generated each time an Atom needs to be
      accessed, as this reduces the number of depend objects that need to be
      held at any one time.

      Args:
         index: The atom that needs to be accessed.
         value: The atom object that holds the new value.
      """

      pat = Atom(self,index)
      pat.p = value.p
      pat.q = value.q 
      pat.m = value.m
      pat.name = value.name

   def get_msum(self):
      """Calculates the total mass."""

      return self.m.sum()
   
   def mtom3(self):
      """Returns a 3*n mass array.

      Returns:
         An array of 3*n elements where each element of m has been copied 
         three times. Used when each degree of freedom needs to be divided 
         by the mass.
      """

      m3 = np.zeros(3*self.natoms,float)
      m3[0:3*self.natoms:3] = self.m
      m3[1:3*self.natoms:3] = m3[0:3*self.natoms:3]
      m3[2:3*self.natoms:3] = m3[0:3*self.natoms:3]
      return m3
                     
   def get_kin(self):
      """Calculates the total kinetic energy of the system."""

      p = depstrip(self.p)
      return 0.5*np.dot(p,p/depstrip(self.m3))
      
   def get_kstress(self):
      """Calculates the total contribution of the atoms to the kinetic stress 
      tensor.
      """

      ks = numpy.zeros((3,3),float)
      ks[0,0] = np.dot(self.px,self.px/self.m)
      ks[1,1] = np.dot(self.py,self.py/self.m)
      ks[2,2] = np.dot(self.pz,self.pz/self.m)
      ks[0,1] = np.dot(self.px,self.py/self.m)
      ks[0,2] = np.dot(self.px,self.pz/self.m)
      ks[1,2] = np.dot(self.py,self.pz/self.m)                        
      return ks
      
class RestartAtoms(Restart):
   """Atoms restart class.

   Handles generating the appropriate atoms class from the xml input file,
   and generating the xml checkpoint tags and data from an instance of the 
   object.

   Attributes:
      natoms: An optional integer giving the number of atoms. Defaults to 0.
      q: An optional array giving the atom positions. Defaults to an empty
         array with no elements.
      p: An optional array giving the atom momenta. Defaults to an empty
         array with no elements.
      m: An optional array giving the atom masses. Defaults to an empty
         array with no elements.
      names: An optional array giving the atom names. Defaults to an empty
         array with no elements.
      from_file: An optional string giving a pdb format file with the atom
         positions. Defaults to ''.
      init_temp: An optional float giving the kinetic temperature to 
         initialise the atom momenta to.
   """

   fields={ "natoms" : (RestartValue, (int,0)), "q" : (RestartArray,(float,np.zeros(0))),  "p" : (RestartArray,(float,np.zeros(0))),
            "m" : (RestartArray,(float, np.zeros(0))),  "names" : (RestartArray,(str,np.zeros(0, np.dtype('|S6')))),
            "from_file" : (RestartValue,(str, "")), "init_temp": (RestartValue, (float, -1.0))  }
       
   def __init__(self, atoms=None, filename=""):
      """Initialises RestartAtoms.

      Args:
         atoms: An optional Atoms object from which to initialise from.
         filename: An optional string giving a filename to take the atom 
            positions from. Defaults to ''.
      """

      super(RestartAtoms,self).__init__()
      if not atoms is None:
         self.store(atoms, filename="")
                       
   def store(self, atoms, filename=""):
      """Takes an Atoms instance and stores a minimal representation of it.

      Args:
         atoms: An Atoms object from which to initialise from.
         filename: An optional string giving a filename to take the atom 
            positions from. Defaults to ''.
      """

      self.natoms.store(atoms.natoms)
      self.q.store(depstrip(atoms.q))
      self.p.store(depstrip(atoms.p))
      self.m.store(depstrip(atoms.m))
      self.names.store(depstrip(atoms.names))
      self.from_file.store(filename)
      
   def fetch(self):
      """Creates an atoms object.

      Returns:
         An atoms object of the appropriate type and with the appropriate
         properties given the attributes of the RestartAtoms object.
      """

      self.check()
      atoms=Atoms(self.natoms.fetch())
      atoms.q=self.q.fetch()
      atoms.p=self.p.fetch()      
      atoms.m=self.m.fetch()   
      atoms.names=self.names.fetch()
      return atoms
   
   def write(self,  name="", indent=""):
      """Overloads Restart write() function so that nothing is written if
      no atoms are present.
      """

      if self.natoms.fetch()>0:
         return super(RestartAtoms,self).write(name=name,indent=indent)
      else:
         return ""
      
   
   def check(self): 
      """Function that deals with optional arguments.

      Deals with the init_temp and from_file arguments, and uses them to 
      intialise some of the atoms parameters depending on which ones have
      been specified explicitly.
      """

      if self.from_file.fetch() != "":
         myatoms, mycell = utils.io.io_pdb.read_pdb(open(self.from_file.fetch(),"r"))
         self.store(myatoms, self.from_file.fetch())      
