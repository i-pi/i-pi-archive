import sys
from engine import *
from utils.io.io_xml import *

def main():
   ifile = open(sys.argv[1],"r")
   xmlrestart = xml_parse_file(ifile)
   
   simrestart = simulation.RestartSimulation()
   simrestart.parse(xmlrestart.fields["simulation"])
   simul = simrestart.fetch()
   
   simul.run()
   del simul

if __name__ == '__main__':
   main()
