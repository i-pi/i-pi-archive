import math, numpy
import atoms
from utils import depend

class Necklace(atoms.Atom):

   def __init__(self, beads, nbeads = 8, temp = 1.0)
      self.beads = beads
      self.temp = depend(name='temp', value=temp)
      self.spring_pot = depend(name='pot', func = self.get_pot, deplist = [self.temp])
      
   def __str__(self):
      rstr = "ATOMS (", str(nbeads), "):\n"
      for i in range(self.nbeads):
         rstr += str(beads[i])
      return rstr

   def pot(self):
      pot = 0.0
      for i in range(self.nbeads-1):
         pot += self.beads[i].pot()
         
         spring_const = self.beads[i].mass/(self.betan*self.hbar)**2
         for j in range(3):
            pot += spring_const*(self.beads[i+1].q[j] - self.beads[i].q[j])**2
         
      pot += self.beads[self.nbeads-1].pot()
      spring_const = self.beads[self.nbeads-1].mass/(self.betan*self.hbar)**2
      for j in range(3):
         pot += spring_const*(self.beads[0].q[j] - self.beads[self.nbeads-1].q[j])**2
     return pot 
   
   def kinetic(self):
      ke = 0.0
      for i in range(self.nbeads):
         ke += self.beads[i].kinetic()
      return ke
