import math, numpy
import atoms

class Necklace(atoms.Atom):

   def __init__(self, qpslice, temp = 1.0, nbeads = 8, name="X", mass=1.0):
      self.nbeads = nbeads
      self.k_Boltz = 1.0
      self.h_bar = 1.0
      self.temp = temp
      self.betan = 1.0/(self.k_Boltz*temp*nbeads)
      self.beads = []
      for i in range(nbeads):
         #this currently creates a reference to the __qp matrix for ALL the different beads. This has to be changed somehow...
         self.beads.append(atoms.Atom(qpslice, name, mass))
      
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
