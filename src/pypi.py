from numpy import *
from engine import *
import sys
from engine import test_Thermo
from engine import io_system

print "hello world"

nat = 3
allthing = zeros((nat,6), float)

f = open("./testfile.txt","r")
syst=engine.System(f)
print syst

syst.step(1.0)
print syst
syst.apply_pbc()
print syst
io_system.print_pdb(syst.atoms, syst.cell)

print syst.kinetic()

#################################

x=allthing[0,0:3]
p=allthing[0,3:6]

print allthing

x[0]=1
p[2]=2
print allthing

x11=allthing[1,1]
x11=4
print allthing

x11=allthing[1:2,1]
x11[0]=4
print allthing

print
print "first cell = ", syst.cell
a, b, c, alpha, beta, gamma = cell.h2abc(syst.cell.h)

print "cell in new coordinates: ", a, b, c, alpha, beta, gamma
syst.cell.h=cell.abc2h(a, b, c, alpha, beta, gamma)

print "back to the start?", syst.cell
print

io_system.print_pdb(syst.atoms,syst.cell)

f.close()
f = open("./testfile.txt", "r")

alist,cell, natoms = io_system.read_pdb(f)

print alist
print cell
print "natoms = ", natoms

myih=syst.cell.ih
myih=syst.cell.ih
myih=syst.cell.ih

print "Trying to call the setter"
print syst.cell.h
print syst.cell.pot()
print syst.cell.pot()

hh = 2*identity(3, float)
syst.cell.h = hh
print syst.cell.pot()
hh[0,1] = -0.2
hh[0,2] = 0.7
hh[1,2] = 1.5
syst.cell.h = hh
print syst.cell.pot()
print "Setter called?"
myih=syst.cell.ih


print syst.cell.h

print syst.atoms[0]
print syst.atoms[1]
print syst.cell.minimum_distance(syst.atoms[0], syst.atoms[1])


print  "before", syst.atoms[0]
print syst.cell
syst.atoms[0] = syst.cell.apply_pbc(syst.atoms[0])
print "after", syst.atoms[0]

lang=langevin.Thermo_Langevin()
lang.dt=4
print lang.dt

print lang.temp

test = test_Thermo.Thermo_test()
test.dt = 4
print test.dt

print test.temp

print "goodbye world"
#print sys.atoms[3].pos.x, sys2.atoms[3].pos.x


