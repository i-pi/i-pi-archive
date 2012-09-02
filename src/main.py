"""Main script from which the simulation is run.

Deals with creation of the simulation object, reading the input file and
initialising the system.

Run using:
      python main.py input_file.xml

Where 'input_file.xml' should be replaced by the name of the xml input_file from
which the system data will be read. For a description of the input file, see the
reference manual.

Functions:
   main: Runs the simulation.
"""

__all__ = ['main']

import sys
from engine import simulation
from inputs.simulation import InputSimulation
from utils.io.io_xml import *

def main(file_name):
   """Runs the simulation.

   Will run automatically when the module is used as a script.
   """

   ifile = open(file_name,"r")
   xmlrestart = xml_parse_file(ifile) # Parses the file.

   simrestart = InputSimulation()
   # Checks the input and partitions it appropriately.
   simrestart.parse(xmlrestart.fields[0][1])

   simul = simrestart.fetch() # Creates the appropriate simulation object.
   simul.run()
   del simul

if __name__ == '__main__':
   main(sys.argv[1])
