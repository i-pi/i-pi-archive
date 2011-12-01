import numpy, math
from barostat import *
from utils.depend import *
from utils import units
import upper_T

class Bussi_S(barostat):
   def __init__(self, pext=numpy.zeros((3,3)), dt = 1.0, w = 1.0, temp = 1.0):
      super(Bussi_S, self).__init__(pext = pext, dt = dt)
      
      self.w = depend(value=w, name='w')
      self.temp = depend(value=temp, name='temp')

   def bind(self, syst):
      self.syst = syst
      self.cell = syst.cell
      syst.cell.w.set(self.w.get())
      syst.cell.pext.set(self.pext.get_array())

      self.kstress = depend(value=numpy.zeros((3,3)), name='kstress', func=self.get_kstress)
      self.piext = depend(value=numpy.zeros((3,3)), name='piext', func=self.get_piext, deplist=[self.pext, self.cell.V, self.cell.V0, self.cell.ih0, self.cell.h])
      self.vir = depend(value=numpy.zeros((3,3)), name='piext', func=self.get_vir, deplist=[self.syst.vir])
      self.stress = depend(value=numpy.zeros((3,3)), name='stress', func=self.get_stress, deplist=[self.vir, self.kstress])

      for atom in self.syst.atoms:
         atom.kstress.add_dependant(self.kstress)

      self.stress.taint(taintme=True); self.piext.taint(taintme=True); self.kstress.taint(taintme=True); self.vir.taint(taintme=True)

   def get_stress(self):
      return self.kstress.get() + self.vir.get()

   def get_kstress(self):
      ks = numpy.zeros((3,3))
      for atom in self.syst.atoms:
         ks += atom.kstress.get()
      ks /= self.cell.V.get()
      return ks

   def get_vir(self):
      return self.syst.vir.get()

   def get_piext(self):
      root = numpy.dot(self.cell.h.get_array(), self.cell.ih0.get_array())
      pi = numpy.dot(root, self.pext.get_array())
      pi = numpy.dot(pi, numpy.transpose(root))
      pi *= self.cell.V0.get()/self.cell.V.get()
      return pi

   def exp_p(self):
      dist_mat = self.cell.p.get_array()*self.dt.get()/self.w.get()
      eig, eigvals = upper_T.compute_eigp(dist_mat)
      i_eig = upper_T.compute_ih(eig)

      exp_mat = numpy.zeros((3,3))
      neg_exp_mat = numpy.zeros((3,3))
      for i in range(3):
         exp_mat[i,i] = math.exp(eigvals[i])
         neg_exp_mat[i,i] = math.exp(-eigvals[i])

      exp_mat = numpy.dot(eig, exp_mat)
      exp_mat = numpy.dot(exp_mat, i_eig)
         
      neg_exp_mat = numpy.dot(eig, neg_exp_mat)
      neg_exp_mat = numpy.dot(neg_exp_mat, i_eig)

      return exp_mat, neg_exp_mat

   def pstep(self):
      V = self.cell.V.get(); dthalf = self.dt.get()/2.0
      L = numpy.zeros((3,3))
      for i in range(3):
         L[i,i] = 3.0 - i

      self.cell.p.get_array()[:] += dthalf*(V*(self.stress.get() - self.piext.get()) + 2.0*units.kb*self.temp.get()*L)

      for i in range(len(self.syst.atoms)):
         atom_i = self.syst.atoms[i]
         self.cell.p.get_array()[:] += dthalf**2/(2.0*atom_i.mass.get())*(numpy.outer(atom_i.f.get_array(), atom_i.p.get_array()) + numpy.outer(atom_i.p.get_array(), atom_i.f.get_array()))
         self.cell.p.get_array()[:] += dthalf**3/(3.0*atom_i.mass.get())*numpy.outer(atom_i.f.get_array(), atom_i.f.get_array())
         self.cell.p.taint(taintme=False)

      self.syst.p.get_array()[:] += self.syst.f.get_array()[:]*dthalf
      self.syst.p.taint(taintme=False)

   def rstep(self):
      """Takes the atom positions, velocities and forces and integrates the 
         equations of motion forward by a step dt"""

      vel_mat = self.cell.p.get_array()/self.w.get()
#      exp_mat, neg_exp_mat = upper_T.Crank_Nicolson(vel_mat*self.dt.get())
      exp_mat, neg_exp_mat = self.exp_p()
      sinh_mat = 0.5*(exp_mat - neg_exp_mat)
      ip_mat = cell.ut_inverse(vel_mat)

      for i in range(len(self.syst.atoms)):
         atom_i = self.syst.atoms[i]
         (atom_i.q.get_array())[:] = numpy.dot(exp_mat, atom_i.q.get_array()) + numpy.dot(ip_mat, numpy.dot(sinh_mat, atom_i.p.get_array()/atom_i.mass.get()))
         (atom_i.p.get_array())[:] = numpy.dot(neg_exp_mat, atom_i.p.get_array())

      self.syst.q.taint(taintme=False);   self.syst.p.taint(taintme=False)          
      (self.cell.h.get_array())[:] = numpy.dot(exp_mat, self.cell.h.get_array())
      self.cell.h.taint(taintme=False)

   def step(self):
      self.pstep()
      self.rstep()
      self.pstep()
