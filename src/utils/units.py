class Constants:
   kb = 3.1668152e-06
   hbar = 1.0
   amu = 1822.8885

class Elements(dict):
   mass_list={
     "X"   :    1.0000/Constants.amu, 
     "Si"  :   28.0860,
     "H"   :    1.0080,
     "H2"  :    2.0160,
     "Li"  :    6.9410,     
     "Ar"   :  39.9480
   }   
   
   @classmethod
   def mass(cls, label):
      return cls.mass_list[label]*Constants.amu

