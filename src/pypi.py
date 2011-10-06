from numpy import *
from engine import *
import sys



print "hello world"

nat = 3
allthing = zeros((nat,6), float)

syst=engine.System(4)
print syst

syst.step(1.0)
print syst

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

io_system.print_pdb(syst.atoms,syst.cell)

alist,cell = io_system.read_pdb(sys.stdin)

myih=syst.cell.ih
myih=syst.cell.ih
myih=syst.cell.ih

print "Trying to call the setter"
hh = 2*identity(3, float)
syst.cell.h = hh
print "Setter called?"
myih=syst.cell.ih

print syst.cell.h

print  "before", syst.atoms[0]
syst.cell.apply_pbc(syst.atoms[0])
print "after", syst.atoms[0]

print "goodbye world"
#print sys.atoms[3].pos.x, sys2.atoms[3].pos.x


