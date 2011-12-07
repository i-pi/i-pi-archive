class Constants:
   kb = 1.0
   hbar = 1.0
   amu = 1822.8885

class Elements(dict):
   mass_list={
     "X"   :    1.0000/Constants.amu, 
     "Si"  :   28.0860,
     "H"   :    1.0080,
     "Ar"   :  39.9480
   }   
   
   def mass(self, label):
      return mass_list[label]*Constants.amu

