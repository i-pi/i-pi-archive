import math, time
import numpy as np
from utils.depend import *
from utils.restart import *
from utils.units import *
from utils.mathtools import eigensystem_ut3x3, invert_ut3x3, exp_ut3x3, det_ut3x3
from engine.thermostats import Thermostat, RestartThermo

class Barostat(dobject): 
   def __init__(self, pext=0.0, sext=None, dt=None, temp=None, thermostat=None):
     
      # This kind of stretches the concept of synced dependencies: sext holds more information than pext but....
      sync_ext=synchronizer()
      dset(self,"sext",depend_array(name='sext', value=np.zeros((3,3)), synchro=sync_ext, func={"pext" : self.p2s} ) )
      dset(self,"pext",depend_value(name='pext', value=0.0, synchro=sync_ext, func={"sext" : self.s2p} ) )            
      if sext is None:
         self.pext=pext
      else:
         self.sext=sext
      

      if thermostat is None:
         thermostat=Thermostat()
      self.thermostat=thermostat   
      # binds options for dt and temperature of the thermostat to those in the barostat
     
      dset(self,"dt",depend_value(name='dt'))
      if dt is None:
         self.dt=2.0*self.thermostat.dt
      else:
         self.dt=dt
      dset(self.thermostat,"dt",   #this involves A LOT of piping
         depend_value(name="dt", func=self.get_halfdt,dependencies=[dget(self,"dt")],dependants=dget(self.thermostat,"dt")._dependants)  )
           
      #temp
      dset(self, "temp", depend_value(name="temp", value=temp))
      deppipe(self, "temp", self.thermostat,"temp")
      if not temp is None: self.temp=temp
      
      self.timer=0.0
      
   def get_halfdt(self):  return self.dt*0.5
         
   def bind(self, beads, cell, forces):
      """Binding function which prepares most of the stuff which will be necessary for a barostat"""
      self.beads=beads
      self.cell=cell
      self.forces=forces

      dset(self,"pot",depend_value(name='pot', func=self.get_pot, 
          dependencies=[ dget(cell,"V0"), dget(cell,"strain"), dget(self,"sext")  ]  ) )            
      dset(self,"piext",depend_value(name='piext', func=self.get_piext, 
          dependencies=[ dget(cell,"V0"), dget(cell,"V"), dget(cell,"h"), dget(cell,"ih0"), dget(cell,"strain"), dget(self,"sext")  ] ) )     
      dset(self,"stress",depend_value(name='stress', func=self.get_stress, 
          dependencies=[ dget(beads,"kstress"), dget(cell,"V"), dget(forces,"vir")  ]  ) )
      dset(self,"press",depend_value(name='press', func=self.get_press, 
          dependencies=[ dget(self,"stress") ] ) )
                
   def s2p(self): return np.trace(self.sext)/3.0
   def p2s(self): return self.pext*np.identity(3)
      
   def pstep(self): pass
   def qstep(self): pass   
      
   def step(self):
      """Dummy atoms barostat step""" 
      self.thermostat.step()
      self.pstep()
      self.qcstep()
      self.pstep()
      self.thermostat.step()
      
   def get_pot(self):
      """Calculates the elastic strain energy of the cell"""
      return self.cell.V0*np.trace(np.dot(self.sext, self.cell.strain))

   def get_piext(self):
      """Calculates the external stress tensor"""
      root = np.dot(depstrip(self.cell.h), depstrip(self.cell.ih0))
      pi = np.dot(root, depstrip(self.sext))
      
      pi = np.dot(pi, np.transpose(root))
      pi *= self.cell.V0/self.cell.V
      return pi
      
   def get_stress(self):
      """Calculates the elastic strain energy of the cell"""
      
      #return (self.beads.kstress+self.forces.vir/self.beads.nbeads)/self.cell.V
      return (np.identity(3)*self.beads.natoms*Constants.kb*self.temp/self.beads.nbeads + self.forces.vir/self.beads.nbeads)/self.cell.V

#TODO  make this something that isn't utter rubbish
#TODO  also include a possible explicit dependence of U on h

   def get_press(self):
      return np.trace(self.stress)/3.0


class BaroFlexi(Barostat):

   def bind(self, beads, cell, forces):
      super(BaroFlexi,self).bind(beads, cell, forces)
      self.thermostat.bind(cell=self.cell)

   def pstep(self):
      
      dthalf = self.dt*0.5
      dthalf2 = dthalf**2/2.0
      dthalf3 = dthalf**3/3.0     

      L = np.zeros((3,3))
      for i in range(3):
         L[i,i] = 3.0 - i
      
      self.cell.p += dthalf*(self.cell.V*(self.stress - self.piext) + 2.0*Constants.kb*self.thermostat.temp*L)       

      m = depstrip(self.beads.m)

      fc = depstrip(self.forces.fnm)/math.sqrt(self.beads.nbeads)
      fx = fc[0:3*self.beads.natoms:3]
      fy = fc[1:3*self.beads.natoms:3]
      fz = fc[2:3*self.beads.natoms:3]
      fxm = fx/m
      fym = fy/m
      fzm = fz/m

      pc = depstrip(self.beads.pc)
      px = pc[0:3*self.atoms.natoms:3]
      py = pc[1:3*self.atoms.natoms:3]
      pz = pc[2:3*self.atoms.natoms:3]
      
      cp = np.zeros((3,3),float)
      cp[0,0] = dthalf2*2.0*np.dot(fxm,px) + dthalf3*np.dot(fx,fxm)
      cp[1,1] = dthalf2*2.0*np.dot(fym,py) + dthalf3*np.dot(fy,fym)
      cp[2,2] = dthalf2*2.0*np.dot(fzm,pz) + dthalf3*np.dot(fz,fzm)
      cp[0,1] = dthalf2*(np.dot(fxm,py) + np.dot(px,fym)) + dthalf3*np.dot(fx,fym)
      cp[0,2] = dthalf2*(np.dot(fxm,pz) + np.dot(px,fzm)) + dthalf3*np.dot(fx,fzm)
      cp[1,2] = dthalf2*(np.dot(fym,pz) + np.dot(py,fzm)) + dthalf3*np.dot(fy,fzm)            
      self.cell.p += cp
      self.beads.p += self.forces.f*dthalf      
      
   def qcstep(self):
      """Takes the atom positions, velocities and forces and integrates the 
         equations of motion forward by a step dt"""
      vel_mat = depstrip(self.cell.p)/self.cell.m

      dist_mat = vel_mat*self.dt
      exp_mat = exp_ut3x3(dist_mat)
      neg_exp_mat = invert_ut3x3(exp_mat)
      sinh_mat = 0.5*(exp_mat - neg_exp_mat)
      ips_mat = np.dot( sinh_mat, invert_ut3x3(vel_mat) )

      pc = depstrip(self.beads.pc).reshape((self.beads.natoms,3)) 
      qc = depstrip(self.beads.qc).reshape((self.beads.natoms,3))
      m3 = depstrip(self.beads.m3)[0,:].reshape((self.beads.natoms,3))       

      qc = np.dot(qc,exp_mat.T)+np.dot(pc/m3,ips_mat.T)
      pc = np.dot(pc,neg_exp_mat.T)

      self.beads.qnm[0,:] = qc.reshape(3*nat)*math.sqrt(self.beads.nbeads)
      self.beads.pnm[0,:] = pc.reshape(3*nat)*math.sqrt(self.beads.nbeads)
                    
      self.cell.h = np.dot(exp_mat, self.cell.h)

      
class BaroRigid(Barostat):

   def get_pot(self):
      """Calculates the elastic strain energy of the cell"""
      return self.cell.V*self.pext
      
   def bind(self, beads, cell, forces):
      super(BaroRigid,self).bind(beads, cell, forces)
      self.thermostat.bind(pm=(self.cell.P, self.cell.M))
      dset(self,"pot",depend_value(name='pot', func=self.get_pot, 
          dependencies=[ dget(self.cell,"V"), dget(self,"pext")  ] ) )
      
   def pstep(self):
      
      dthalf = self.dt*0.5
      dthalf2 = dthalf**2/2.0
      dthalf3 = dthalf**3/3.0     
      
      self.cell.P += dthalf*3.0*(self.cell.V*(self.press - self.pext) + 2.0*Constants.kb*self.temp)

      fc = depstrip(self.forces.fnm)[0,:]/math.sqrt(self.beads.nbeads)
      m = depstrip(self.beads.m3)[0,:]
      pc = depstrip(self.beads.pc)
            
      self.cell.P += dthalf2*np.dot(pc,fc/m) + dthalf3*np.dot(fc,fc/m)
   
      self.beads.p += depstrip(self.forces.f)*dthalf      
           
   def qcstep(self):
      """Takes the atom positions, velocities and forces and integrates the 
         equations of motion forward by a step dt"""
      vel = self.cell.P[0]/self.cell.m
      exp, neg_exp = (math.exp(vel*self.dt), math.exp(-vel*self.dt))
      sinh = 0.5*(exp - neg_exp)

      pc = depstrip(self.beads.pc)
      qc = depstrip(self.beads.qc)
      m = depstrip(self.beads.m3)[0,:]      
      qc*=exp
      qc+=(sinh/vel)* pc/m
      pc *= neg_exp

      self.beads.qnm[0,:] = qc*math.sqrt(self.beads.nbeads)
      self.beads.pnm[0,:] = pc*math.sqrt(self.beads.nbeads)

      self.cell.V*=exp**3

      
class RestartBaro(Restart):
   attribs={ "kind": (RestartValue, (str, "rigid")) }
   fields={ "thermostat": (RestartThermo, ()) }
   
   def store(self, baro):
      if type(baro) is BaroRigid:
         self.kind.store("rigid")
      if type(baro) is BaroFlexi:
         self.kind.store("flexible")
      else:
         self.kind.store("unknown")      
      self.thermostat.store(baro.thermostat)
      
   def fetch(self):
      if self.kind.fetch().upper() == "RIGID":
         baro=BaroRigid(thermostat=self.thermostat.fetch())
      elif self.kind.fetch().upper() == "FLEXIBLE":
         baro=BaroFlexi(thermostat=self.thermostat.fetch())
      else:
         baro=Barostat(thermostat=self.thermostat.fetch())

      return baro
