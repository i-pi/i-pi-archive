from numpy import *
from engine import *
import sys
from engine import io_system
from engine import dynamics
from engine import forces
from engine import rp_engine
from engine import rp_dynamics
from engine import PILE
from engine import Bussi
import numpy, random
from utils.mathtools import *
from utils.mass_list import *

print "hello world"

#a = numpy.zeros((3,3),float)
#b = numpy.array(range(9)); b.shape=(3,3)

#c=a[:]
#c[:]=b
#a[1,1]=-1
#print "a", a
#print "b", b

#exit()

#ll=zeros((12,22),float)
#aa=engine.Atom(ll)
#aa.q[0]=2
#print aa.q[0], aa.q._dep_proxy__dep._depend__tainted
#exit()

#nat = 3
#allthing = zeros((nat,6), float)
#f = open("./testfile.txt","r")
#syst=engine.System.from_pdbfile(f)


##print "Get atomic ke"
##for i in range(0,syst.natoms):
##   print "KE(",i,"): ",syst.atoms[i].kin
##print "Global KE:", syst.kin
##
##print "Now we modify the momentum of atom 0, and get again ke's"
##syst.atoms[0].p[0]=10
##for i in range(0,syst.natoms):
##   print "KE(",i,"): ",syst.atoms[i].kin
##print "Global KE:", syst.kin
##
##print "Now we modify the global momentum, and get again ke's"
##syst.p[2]=1
##for i in range(0,syst.natoms):
##   print "KE(",i,"): ",syst.atoms[i].kin
##print "Global KE:", syst.kin
##
##print "Now we modify the momentum of atom 0, and 3 and get just the global ke"
##syst.atoms[0].p[1]=10
##syst.atoms[3].p[2]=10
##print "Global KE:", syst.kin

##exit()
##syst.step(1.0)
##print syst
##syst.apply_pbc()
##print syst
##io_system.print_pdb(syst.atoms, syst.cell)

##print syst.kinetic()

##################################

#x=allthing[0,0:3]
#p=allthing[0,3:6]

#print allthing

#x[0]=1
#p[2]=2
#print allthing

#x11=allthing[1,1]
#x11=4
#print allthing

#x11=allthing[1:2,1]
#x11[0]=4
#print allthing

#print
#print "first cell = ", syst.cell
#a, b, c, alpha, beta, gamma = cell.h2abc(syst.cell.h)

#print "cell in new coordinates: ", a, b, c, alpha, beta, gamma
#syst.cell.h=cell.abc2h(a, b, c, alpha, beta, gamma)

#print "back to the start?", syst.cell
#print

#io_system.print_pdb(syst.atoms,syst.cell)

#f.close()
#f = open("./testfile.txt", "r")

#alist,cell, natoms = io_system.read_pdb(f)

#print alist
#print cell
#print "natoms = ", natoms

#myih=syst.cell.ih
#myih=syst.cell.ih
#myih=syst.cell.ih

#print "Trying to call the setter"
#print syst.cell.h
#print syst.cell.pot()
#print syst.cell.pot()

#hh = 2*identity(3, float)
#syst.cell.h = hh
#print syst.cell.pot()
#hh[0,1] = -0.2
#hh[0,2] = 0.7
#hh[1,2] = 1.5
#syst.cell.h = hh
#print syst.cell.pot()
#print "Setter called?"
#myih=syst.cell.ih


#print syst.cell.h

#print syst.atoms[0]
#print syst.atoms[1]
#print syst.cell.minimum_distance(syst.atoms[0], syst.atoms[1])


#print  "before", syst.atoms[0]
#print syst.cell
#syst.atoms[0] = syst.cell.apply_pbc(syst.atoms[0])
#print "after", syst.atoms[0]

#lang=langevin.Thermo_Langevin()
#lang.dt=4
#print lang.dt

#print lang.temp

#test = test_Thermo.Thermo_test()
#test.dt = 4
#print test.dt

#print test.temp
#print


#f = "engine/system_out.xml"
#ffield = forces.forcefield()
#g = open("./testfile2.txt", "r")
#syst = engine.System.from_pdbfile(g)
#ffield.bind(syst.cell, syst.atoms, syst.pot, syst.f, syst.vir)
#
#io_system.xml_read(f, ffield)
#
#print syst.f.get()
#print
#print ffield.f.get()
#print
#print syst.vir.get()
#print
#print ffield.vir.get()
#print
#print syst.pot.get()
#print
#print ffield.pot.get()
#print
#
#g.close()

f = open("./fcp4cell.pdb", "r")
#f = open("./testfile.txt", "r")
#syst = engine.System.from_pdbfile(f)
#g = open("./forces/system.xml", "w")
#io_system.xml(syst, g)
#print syst
#exit()
#thermo = langevin.Thermo_Langevin(dt = 0.1)

#syst=rp_engine.RP_sys.from_pdbfile(f, forces.rp_pipeforce( {"pipein": "forces/pipeforce", "pipeout": "forces/pipepos"}), nbeads = 12, temp = 25/3.1577464/10**5 )
#thermo = langevin.langevin(tau=1e-1)
#thermo = PILE.PILE(tau_0 = 1e3)
#for system in syst.systems:
#   system.cell.w.set(1e1)
#
#syst.cell.w.set(1e1)
#
##nvt = rp_dynamics.rp_nvt_ensemble(syst=syst, thermo=thermo, dt = 206.706895, temp=25/3.1577464/10**5)
#nvt = rp_dynamics.rp_nvt_ensemble(syst=syst, thermo=thermo, dt = 0.100*206.706895, temp=25/3.1577464/10**5)
##nvt = dynamics.nve_ensemble(syst=syst, dt = 2.5e-4)
#print "#Initial vir is ", syst.vir.get()
##print "#Initial f.at is ", syst.systems[0].atoms[0].f.get_array()
##print "#Initial f is ", syst.f.get_array()[0,0:3]
##print "#Final f.at is ", syst.systems[0].atoms[0].f.get_array()
##print "#Initial cell p is ", syst.cell.p.get()
#print "# Initial pot is ", syst.pot.get()
#print "# Initial ke is ", syst.kin.get(), syst.kin_estimator.get(), 1.5*len(syst.atoms)/(syst.betan.get()*len(syst.systems))
#print "# Initial econs is ", nvt.econs.get()
##print "# Thermo T is ", nvt.thermo.T.get()
#print "# V K ECNS V"
##f = open("./traj6.pdb", "w")
#for istep in range(90):
##for istep in range(1800*10):
#   nvt.step()
#   #io_system.print_pdb_RP(syst.systems, f)
#   print syst.pot_estimator.get(), syst.kin_estimator.get(), nvt.econs.get()
#   q_tilde = numpy.dot(syst.trans_mat.get_array(), syst.q.get_array())
#   
#
#exit()
#print "Equilibration done: starting actual run..."
#kin = 0.0
#for istep in range(2000):
#   nvt.step()
#   
#   io_system.print_pdb_RP(syst.systems, f)
#
##   print syst.pot.get(), syst.kin.get(), nvt.econs.get(), syst.systems[0].cell.V.get()
#   print syst.pot_estimator.get(), syst.kin_estimator.get(), nvt.econs.get(), syst.systems[0].cell.V.get()
#   kin += syst.kin_estimator.get()
#kin /= 2000
#print kin
#
#h = open("./forces/pipepos","w")
#io_system.xml_terminate(h)
#h.close()
#
#exit()

pext = numpy.zeros((3,3))
w = mlist.masses["  Ar"] * 256
temp = 1.1761e-4
dt = 445.0
syst=engine.System.from_pdbfile(f, forces.pipeforce( {"pipein": "forces/pipeforce", "pipeout": "forces/pipepos"} ), w = w, pext = pext )
#syst=engine.System.from_pdbfile(f, forces.LJ( {"eps": 0.1, "sigma": 0.38, "rc": 0.38*2.5} ) )
thermo = langevin.langevin(temp = temp, dt = dt/2, tau=1e3)
thermo_cell = langevin.langevin(temp = temp, dt = dt/2, tau=1e3)
#syst.cell.w.set(1e1)
#syst.cell.w.set(mlist.masses["  Ar"]*256)
#nvt=dynamics.npt_ensemble(syst=syst, thermo=thermo, cell_thermo=thermo_cell, dt=1e-2, temp=1e-2, pext=10.0)
#pext=10.0*numpy.identity(3); pext[0,2]=pext[2,0]=0
#pext = numpy.zeros((3,3),float)
#nvt=dynamics.nst_ensemble(syst=syst, thermo=thermo, cell_thermo=thermo_cell, dt=445, temp=1.1761e-4)
nvt = dynamics.test_NST(syst = syst, barostat = Bussi.Bussi_S(pext = pext, dt = dt/2.0, w = w, temp = temp), thermo = thermo, cell_thermo = thermo_cell, temp = temp, dt = dt)
#nvt=dynamics.nvt_ensemble(syst=syst, thermo=thermo, dt=100, temp=1.1761e-4)

print "#Initial vir is ", syst.vir.get()
#print "#Initial f is ", syst.f.get()
#print "#Initial cell p is ", syst.cell.p.get()
print "# Initial pot is ", syst.pot.get()
print "# Initial econs is ", nvt.econs.get()
print "# Thermo T is ", nvt.thermo.T.get()
print "# V K ECNS V"

c = numpy.zeros((3,3,3,3))
c_bar = numpy.zeros((3,3,3,3))
C = numpy.zeros((6,6))
vol = 0.0
steps = 30000
steps2 = 4000
steps = 400
steps2 = 40

h0_bar = numpy.array(syst.cell.h.get_array())
#f = open("./traj9.pdb", "w")
#g = open("./corr_file2.txt", "w")
for istep in range(steps2):
   nvt.step()
 #  if istep%20 == 0:
  #    io_system.print_pdb(syst.atoms, syst.cell, f)
   print "step ", istep + 1, " of: ", steps + steps2
   print syst.pot.get(), syst.kin.get(), nvt.econs.get(), syst.cell.V.get()
   h0_bar += syst.cell.h.get_array()
   syst.cell.h0.set(h0_bar/(istep+2.0))
#exit()
   
for istep in range(steps):
   nvt.step()
#   if istep%20 == 0:
#      io_system.print_pdb(syst.atoms, syst.cell, f)

   print "step ", steps2 + istep + 1, " of: ", steps + steps2
   print syst.pot.get(), syst.kin.get(), nvt.econs.get(), syst.cell.V.get()
   h0_bar += syst.cell.h.get_array()
   syst.cell.h0.set(h0_bar/(steps2 + istep + 2.0))
   vol += syst.cell.V.get()
   for i in range(3):
      for j in range(3):
         for k in range(3):
            for l in range(3):
               c_bar[i,j,k,l] += syst.cell.strain.get()[i,j]*syst.cell.strain.get()[k,l]
               c[i,j,k,l] = syst.cell.strain.get()[i,j]*syst.cell.strain.get()[k,l]


   C[0,0] = c[0,0,0,0]
   C[1,1] = c[1,1,1,1]
   C[2,2] = c[2,2,2,2]
   C[0,1] = C[1,0] = c[0,0,1,1]
   C[0,2] = C[2,0] = c[0,0,2,2] 
   C[2,1] = C[1,2] = c[2,2,1,1] 
   C[0,3] = C[3,0] = c[0,0,1,2]*2.0
   C[0,4] = C[4,0] = c[0,0,2,0]*2.0
   C[0,5] = C[5,0] = c[0,0,0,1]*2.0
   C[1,3] = C[3,1] = c[1,1,1,2]*2.0
   C[1,4] = C[4,1] = c[1,1,2,0]*2.0
   C[1,5] = C[5,1] = c[1,1,1,0]*2.0
   C[2,3] = C[3,2] = c[2,2,1,2]*2.0
   C[2,4] = C[4,2] = c[0,2,2,2]*2.0
   C[2,5] = C[5,2] = c[2,2,0,1]*2.0
   C[3,3] = c[1,2,1,2]*4.0
   C[4,4] = c[0,2,0,2]*4.0
   C[5,5] = c[1,0,1,0]*4.0
   C[3,4] = C[4,3] = c[1,2,2,0]*4.0
   C[3,5] = C[5,3] = c[1,2,1,0]*4.0
   C[5,4] = C[4,5] = c[1,0,2,0]*4.0

   #g.write(str((C[0,0] + C[1,1] + C[2,2])/(3.0)) + "  "), g.write(str((C[0,1] + C[0,2] + C[1,2])/(3.0)) + "  "), g.write(str((C[4,4] + C[5,5] + C[3,3])/(3.0)) + "\n")

C[0,0] = c_bar[0,0,0,0]
C[1,1] = c_bar[1,1,1,1]
C[2,2] = c_bar[2,2,2,2]
C[0,1] = C[1,0] = c_bar[0,0,1,1]
C[0,2] = C[2,0] = c_bar[0,0,2,2] 
C[2,1] = C[1,2] = c_bar[2,2,1,1] 
C[0,3] = C[3,0] = c_bar[0,0,1,2]*2.0
C[0,4] = C[4,0] = c_bar[0,0,2,0]*2.0
C[0,5] = C[5,0] = c_bar[0,0,0,1]*2.0
C[1,3] = C[3,1] = c_bar[1,1,1,2]*2.0
C[1,4] = C[4,1] = c_bar[1,1,2,0]*2.0
C[1,5] = C[5,1] = c_bar[1,1,1,0]*2.0
C[2,3] = C[3,2] = c_bar[2,2,1,2]*2.0
C[2,4] = C[4,2] = c_bar[0,2,2,2]*2.0
C[2,5] = C[5,2] = c_bar[2,2,0,1]*2.0
C[3,3] = c_bar[1,2,1,2]*4.0
C[4,4] = c_bar[0,2,0,2]*4.0
C[5,5] = c_bar[1,0,1,0]*4.0
C[3,4] = C[4,3] = c_bar[1,2,2,0]*4.0
C[3,5] = C[5,3] = c_bar[1,2,1,0]*4.0
C[5,4] = C[4,5] = c_bar[1,0,2,0]*4.0
   
C /= steps
C_inv = numpy.linalg.inv(C)

#c /= steps
#c = 1.0/c
vol /= steps
C_inv = nvt.thermo.temp.get()*units.kb/vol*C_inv
unit = len(syst.atoms)*units.kb*nvt.thermo.temp.get()/vol
#print (c[0,0,0,0] + c[1,1,1,1] + c[2,2,2,2])/(2.0*unit)
#print (c[0,0,1,1] + c[0,0,2,2] + c[1,1,2,2])/(2.0*unit)
#print (c[1,2,1,2] + c[2,0,2,0] + c[0,1,0,1])/(2.0*unit)
print (C_inv[0,0] + C_inv[1,1] + C_inv[2,2])/(3.0*unit)
print (C_inv[0,1] + C_inv[0,2] + C_inv[1,2])/(3.0*unit)
print (C_inv[4,4] + C_inv[5,5] + C_inv[3,3])/(3.0*unit)
print unit

h = open("./forces/pipepos","w")
io_system.xml_terminate(h)
h.close()

pot_func = forces.LJ
kwargs = {"eps": 0.1, "sigma": 0.3, "rc": 0.3*2.5}


#syst2 = dynamics.NST_ens.from_pdbfile(f, thermo, pot_func, temp=1e-2, tau=1e-1,  dt = 0.005, **kwargs)

#print syst2.syst
#print syst2.thermo.dt
#print syst2.thermo.temp

#print 
#print "SIMULATION STARTS HERE!"

#syst2.simulation(100)


#syst2.simulation(1)

#syst3 = dynamics.NST_ens.from_ensemble(syst2)
#syst3 = engine.System.from_system(syst2.syst)
#print
#print "syst2: ", syst2.syst
#print "syst3: ", syst3.syst

#print
#print syst2.syst.cell.p*syst2.dt/syst2.syst.cell.w
#print syst2.exp_p()
#print

print "goodbye world"
#print sys.atoms[3].pos.x, sys2.atoms[3].pos.x


