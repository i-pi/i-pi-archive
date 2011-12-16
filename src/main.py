from engine import *
from utils.depend import *
from utils.io.io_pdb import *
import numpy, time
import pdb


#atoms, cell = read_pdb(open("testfile4.pdb","r"))
#atoms, cell = read_pdb(open("sysfile2.pdb","r"))
##atoms.p = np.random.randn(atoms.natoms*3)*0.5

myatoms, cell = read_pdb(open("sysfile2.pdb","r"))

ratoms=atoms.RestartAtoms()
ratoms.store(myatoms)

force=forces.FFSocket(" 0.0003793865e0  6.43452e0   46.651d0   50.00  ")
force.bind(myatoms, cell)
rforce=forces.RestartForce()
rforce.store(force)

print rforce.write("force")
exit(1);
      
#force=forces.FFLennardJones({"sigma":0.8, "eps":1.0, "rc":2.5*0.8})

#cell.h0=[[0.5,0.25,1],[0,0.1,1],[0,0,1]]


#ens=ensembles.NSTEnsemble(dt=0.01, temp=0.5, thermostat=thermostats.ThermoLangevin(tau=1.0), 
#                          barostat=barostats.BaroFlexi( pext=1.0, thermostat=thermostats.ThermoLangevin(tau=1.0e0) ),
#                          fixcom=True )

print "cell mass is computed as", cell.m
#ens=ensembles.NSTEnsemble(dt=500, temp=1e-4, thermostat=thermostats.ThermoLangevin(tau=1.0e4), 
#                          barostat=barostats.BaroRigid( pext=1e-8, thermostat=thermostats.ThermoLangevin(tau=1.0e4) ),
#                          fixcom=True ); 
ens=ensembles.NVTEnsemble(dt=1000, temp=5e-4, thermostat=thermostats.ThermoLangevin(tau=1.0e4), fixcom=False )
#ens=ensembles.NVTEnsemble(dt=0.01, temp=0.5, thermostat=thermostats.ThermoLangevin(tau=1.0), fixcom=False )
cell=CellRigid(h=cell.h, m=cell.m) #need a flexible cell, actually
#cell.h0=[[0.5,0.25,1],[0,0.1,1],[0,0,1]]

print "USING ", cell.h0
test=simulation.Simulation(atoms, cell, force, ens)

ftraj=open("trajfile.pdb","w")
fout=open("output.dat","w")
#pdb.set_trace()
start=time.clock()
wall=time.time()
for i in range(1000):
#   try:
      test.ensemble.step()

#      fout.write(str(i)+" "+str(test.ensemble.econs)+" "+str(test.atoms.kin)+" "+str(test.force.pot)+" "+str(test.ensemble.thermostat.ethermo)+" "+str(ens.barostat.pot)+" "+str(cell.kin)+" "+str(cell.V)+" "+str(cell.h[0,0]/cell.h[1,1])+"\n")
      fout.write(str(i)+" "+str(test.ensemble.econs)+" "+str(test.atoms.kin)+" "+str(test.force.pot)+" "+str(test.ensemble.thermostat.ethermo)+" "+"\n")
      if i%20 ==0:    print_pdb(test.atoms, test.cell, ftraj)
#   except:
#      test.force.socket.server.shutdown(2)   
#      test.force.socket.server.close()
#      del test
#      exit(1)
print "Time elapsed:", time.clock()-start, "wall time", time.time()-wall
print "of which:", force.timer , " spent in", force.ncall ," force evaluations", "wall time:", force.twall

#print "Time in barostat select", ens.barostat.timer
print "Time in ensemble select", ens.timer

del test.force.socket

#v1=depend_array(value=numpy.random.randn(5,5),name="v1")
#v2=depend_array(value=numpy.random.randn(5),name="v2")
#print v2
#print "created v1 and v2", v1.name, v2.name
#print "creating a view"
#v3=v2[1:4]
#v4=v3[1:3]

#print "v3 name is ",v3.name, v3
#print "v4 name is ",v4.name, v4, v4[1]

#print "setting a view"
#v4[:]=1.0

#print v2

#v3=v1[1]+1.0

#print "creating a reference"
#v4=numpy.dot(v1,v2-v3)
#print "@@V4", type(v4), v4.name, v4

#print "creating a reference"
#v4=v3*v2
#print "@@V4",  type(v4), v4.name, v4.base

#print "adding to the reference"
#v4+=100.0
#print "@@V4",  type(v4), v4.name, v4.base

#print "creating a reference"
#v2+=v3
#print "@@V2",  type(v2), v2.name

#exit(1)
#print "@@INIT CELL"

