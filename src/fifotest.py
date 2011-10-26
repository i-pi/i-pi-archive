import sys, os
import numpy, random
from engine import engine, forces, langevin, dynamics

f = open("./testfile4.txt", "r")

syst=engine.System.from_pdbfile(f, forces.LJ( {"eps": 0.1, "sigma": 0.866, "rc": 0.866*2.5} ) )
thermo = langevin.langevin(tau=1e-1)
thermo_cell = langevin.langevin(tau=1e-1)
syst.cell.w.set(1e1)
pext=10.0*numpy.identity(3); pext[0,2]=pext[2,0]=10
nst=dynamics.nst_ensemble(syst=syst, thermo=thermo, cell_thermo=thermo_cell, dt=4e-3, temp=1e-2, pext=pext)

#print "#Initial vir is ", syst.vir.get()
#print "#Initial f is ", syst.f.get()
#print "# Initial pot is ", syst.pot.get()
#print "# Thermo T is ", nst.thermo.T.get()
#print "# V K ECNS V"

os.mkfifo('fifo1')

g = open('fifo1', 'r')
f = open('fifo1', 'w')
xml_write(nst.syst, f)
xml_read(g, nst.syst.ffield)
f.close()
g.close()
exit()
