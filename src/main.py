from engine import *
from utils.depend import *
from utils.io.io_pdb import *
from utils.io.io_xml import *

import numpy, time, sys
import pdb


ifile=open(sys.argv[1],"r")

buf="";
for a in ifile: buf+=a
xmlrestart=parse_xml(buf)

simrestart=simulation.RestartSimulation();
simrestart.parse(xmlrestart.fields["simulation"])
simul=simrestart.fetch()

simul.run()

fcheck=open("checkpoint.xml","w")
simrestart.store(simul)
del simul.force.socket

fcheck.write(simrestart.write("simulation"))
