from numpy import *
from engine import *


print "hello world"

nat = 3
allthing = zeros((nat,6), float)

syst=atoms.System(4)
print syst

syst.step(1.0)
print syst

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


print "goodbye world"
#print sys.atoms[3].pos.x, sys2.atoms[3].pos.x


