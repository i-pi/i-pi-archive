"""Contains fundamental constants in atomic units.

Classes:
   Constants: Class whose members are fundamental contants.
   Elements: Class which contains the mass of different elements
   Units: Class which contains the methods needed to transform
      between different systems of units.
"""

__all__ = ['Constants', 'Elements', 'Units', 'UnitMap']

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
         Has the form {"label": Mass in a.m.u.}. Note that the generic "X"
         label is assumed to be an electron.
   """

   mass_list={
     "X"   :    1.0000/Constants.amu, 
     "H"   :    1.0079,
     "D"   :    2.0141,
     "Z"   :  1.382943, #an interpolated H--D atom, based on y=1/sqrt(m) scaling
     "Si"  :   28.0860,
     "O"   :   15.9994,
     "Cl"  :   35.453,
     "Na"  :   22.9898,
     "N"   :   14.0067,
     "H2"  :    2.0160,
     "Li"  :    6.9410,     
     "Ar"  :   39.9480,
     "I"   :  126.9045 
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

UnitMap = { 
   "undefined" : {
      "" : 1.00 
      },      
   "energy" :   {
      "" : 1.00,
      "atomic_unit"  : 1.00, 
      "electronvolt" : 0.036749326,
      "joule"        : 2.2937128e+17,
      "kj_mol"       : 0.00038087989,
      "kcal_mol"     : 0.0015946679,
      "kelvin"       : 3.1668152e-06
      },
   "time" :     {
      "" : 1.00,
      "atomic_unit"  : 1.00,
      "second"       : 4.1341373e+16,
      "picosecond"   : 41341.373,
      "femtosecond"  : 41.341373
      },
   "length" :     {
      "" : 1.00,
      "atomic_unit"  : 1.00,
      "angstrom"     : 1.8897261,
      "nanometer"    : 0.18897261
      },
   "velocity":    {
      "" : 1.00,
      "atomic_unit"  : 1.00
      },           
   "force":       {
      "" : 1.00,
      "atomic_unit"  : 1.00
      },
   "mass":        {
      "" : 1.00,
      "atomic_unit"  : 1.00,
      "electronmass" : 1.00/Constants.amu
      },
   "pressure" :     {
      "" : 1.00,
      "atomic_unit"  : 1.00,
      "bar"          : 3.398827377e-9,
      "atmosphere"   : 3.44386184e-9,
      "pascal"       : 3.398827377e-14,
      "kilopascal"   : 3.398827377e-11
      }
}


class Units(object):
   """A class to perform simple unit conversions based on a map.

   Uses the UnitMap object, which is a dictionary of the form 
   {'unit type': {'name': value}}, where the value is the value of the 
   particular unit described by name in atomic units.

   Attributes:
      family: The unit type. Will be either a base unit such as 'time', a 
         derived unit such as 'force', or unitless, which is given by 'number'. 
   """
   
   def __init__(self, family="number"):
      """Initialises Units.

      Args:
         family: The unit type. Defaults to 'number'.

      Raises:
         IndexError: Raised if the unit type given has not been defined.
      """

      if not (family == "number"  or family in UnitMap):
         raise IndexError(family + " is an undefined units kind.")
      self.family = family
      
   def to_internal(self, number, unit):
      """Converts number from the unit type given to internal units using the 
      appropriate conversion factor.

      Args:
         number: The value of the number in the unit given.
         unit: The particular unit of the appropriate unit type that the 
            number is given in.

      Raises:
         TypeError: Raised if the unit given is not of the correct type.
      """

      if self.family == "number":
         return number
      if not unit in UnitMap[self.family]:
         raise TypeError(unit + " is an undefined unit for kind " + self.family + ".")
      return number*UnitMap[self.family][unit]

   def to_user(self, number, unit):
      """Converts number from internal units to the specified unit.

      Args:
         number: The value of the number in internal units.
         unit: The unit for the number to be converted into.
      """

      return number/UnitMap[self.family][unit]
