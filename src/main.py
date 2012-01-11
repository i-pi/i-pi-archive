from engine import *
from utils.depend import *
from utils.io.io_pdb import *
from utils.io.io_xml import *

import numpy, time, sys
import pdb


ifile=open(sys.argv[1],"r")
xmlrestart=xml_parse_file(ifile)

simrestart=simulation.RestartSimulation();
simrestart.parse(xmlrestart.fields["simulation"])
simul=simrestart.fetch()

import utils.depgraph as udg
udg.plot_deps()

pdb.set_trace()
simul.run()

del simul
