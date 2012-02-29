"""Contains fundamental constants in atomic units.

Classes:
   Constants: Class whose members are fundamental contants.
   Elements: Class which contains the mass of different elements
"""

__all__ = ['Constants', 'Elements']

class Constants:
   """Class whose members are fundamental contants.

   Attributes:
      kb: Boltzmann constant.
      hbar: Reduced Planck's constant.
      amu: Atomic mass unit.
   """

   kb = 3.1668152e-06
   hbar = 1.0
   amu = 1822.8885


class Elements(dict):
   """Class which contains the mass of different elements.

   Attributes:
      mass_list: A dictionary containing the masses of different elements.
         Has the form {"label": Mass in a.m.u.}.
   """

   mass_list={
     "X"   :    1.0000/Constants.amu, 
     "Si"  :   28.0860,
     "H"   :    1.0080,
     "O"   :   15.9994,
     "H2"  :    2.0160,
     "Li"  :    6.9410,     
     "Ar"  :   39.9480
   }   
   
   @classmethod
   def mass(cls, label):
      """Function to access the mass_list attribute.
   
      Note that this does not require an instance of the Elements class to be 
      created, as this is a class method. Therefore using Elements.mass(label) 
      will give the mass of the element with the atomic symbol given by label.

      Args:
         label: The atomic symbol of the atom whose mass is required.

      Returns:
         A float giving the mass of the atom with atomic symbol label.
      """

      return cls.mass_list[label]*Constants.amu

