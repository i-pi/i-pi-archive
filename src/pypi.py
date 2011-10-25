from numpy import *
from engine import *
import sys
from engine import io_system
from engine import dynamics
from engine import forces
import numpy, random
from utils.mathtools import *

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

f = "engine/system_out.xml"
ffield = forces.forcefield()
g = open("./testfile2.txt", "r")
syst = engine.System.from_pdbfile(g)
ffield.bind(syst.cell, syst.atoms, syst.pot, syst.f, syst.vir)

io_system.xml_read(f, ffield)

print syst.f.get()
print
print ffield.f.get()
print
print syst.vir.get()
print
print ffield.vir.get()
print
print syst.pot.get()
print
print ffield.pot.get()
print

g.close()

f = open("./testfile2.txt", "r")
#syst = engine.System.from_pdbfile(f)
#g = open("./forces/system.xml", "w")
#io_system.xml(syst, g)
#print syst
#exit()
#thermo = langevin.Thermo_Langevin(dt = 0.1)

syst=engine.System.from_pdbfile(f, forces.LJ( {"eps": 0.1, "sigma": 0.19, "rc": 0.19*2.5} ) )
thermo = langevin.langevin(tau=1e-1)
thermo_cell = langevin.langevin(tau=1e-1)
syst.cell.w.set(1e1)
#nvt=dynamics.npt_ensemble(syst=syst, thermo=thermo, cell_thermo=thermo_cell, dt=1e-2, temp=1e-2, pext=10.0)
pext=10.0*numpy.identity(3); pext[0,2]=pext[2,0]=10
nvt=dynamics.nst_ensemble(syst=syst, thermo=thermo, cell_thermo=thermo_cell, dt=1e-2, temp=1e-2, pext=pext)

print "#Initial vir is ", syst.vir.get()
print "#Initial f is ", syst.f.get()
print "# Initial pot is ", syst.pot.get()
print "# Thermo T is ", nvt.thermo.T.get()
print "# V K ECNS V"
f = open("./traj.pdb", "w")
for istep in range(500):
   nvt.step()
   io_system.print_pdb(syst.atoms, syst.cell, f)

   nvt.econs.get()

   print syst.pot.get(), syst.kin.get(), nvt.econs.get(), syst.cell.V.get()

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


