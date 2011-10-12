import math, numpy

class LJ:

   def __init__(self, syst, eps = 1.0, sigma = 1.0, rc = 2.5)
      self.syst = syst
      self.eps = eps 
      self.sigma = sigma
      self.rc = rc

   def separation(self, atom_i, atom_j):
      rij = self.syst.cell.minimum_distance(atom_i, atom_j)
      r = math.sqrt(numpy.dot(rij, rij))
      return r, rij

   def LJ_force(self, r):
      if (r > self.rc):
         return 0.0
      else:
         return 4*self.eps*(12/r*(self.sigma/r)**12 - 6/r*(self.sigma/r)**6)

   def LJ_fij(self, atom_i, atom_j):
      fij = numpy.zeros(3, float)
      fji = numpy.zeros(3, float)
      r, rij = self.separation(atom_i, atom_j)
      f_tot = self.LJ_force(r)
      for i in range(3):
         fij[i] = f_tot*rij[i]/r
      fji = -fij
      return fij, fji, r

   def LJ_pot(self, r):
      if (r > self.rc):
         return 0.0
      else:
         return 4*self.eps*((self.sigma/r)**12 - (self.sigma/r)**6 - (self.sigma/self.rc)**12 + (self.sigma/self.rc)**6)

   def syst_update(self):
      V = self.syst.cell.V
      natoms = self.syst.natoms

      self.syst.strain = numpy.zeros((3,3),float)
      self.syst.pot = 0.0
      self.syst.f = numpy.zeros(3*natoms, float)
   
      for i in range(natoms):
         atom_i = self.syst.atoms[i]
         p_i = atom_i.p
         mass_i = atom_i.mass

         self.syst.strain += numpy.outer(p_i, p_i)/(mass_i*V)

         for j in range(i+1, natoms):
            atom_j = self.syst.atoms[j]

            fij, fji, r = self.LJ_fij(atom_i, atom_j)
            atom_i.f += fij
            atom_j.f += fji
            self.syst.pot += self.LJ_pot(rij)

            self.syst.strain += numpy.outer(fij, rij)/V
