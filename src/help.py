import sys
from inputs import *
from utils.io.io_xml import *
from optparse import OptionParser

objects = {'barostats': barostats.InputBaro(), 'cell': cell.InputCell(), 'simulation': simulation.InputSimulation(), 'ensembles': ensembles.InputEnsemble(), 'thermostats': thermostats.InputThermo(), 'interface': interface.InputInterface(), 'forces': forces.InputForce(), 'atoms': atoms.InputAtoms(0), 'beads': beads.InputBeads(0,1), 'prng': prng.InputRandom()}

usage = "usage: python %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option("-x", action="store_true", dest = "xml", default=False, help="write an xml help file")
parser.add_option("-l", action="store_true", dest = "latex", default=False, help="write a latex help file")
parser.add_option("-n", action="store", type="int", dest="levels", help="number of levels depth to which data is printed out")
parser.add_option("-i", action="store", dest="opt", help="Root object for the help files. Options: ['barostats', 'cell', 'simulation', 'ensembles', 'thermostats', 'interface', 'forces', 'atoms', 'beads', 'prng']", default='simulation')
(options, args) = parser.parse_args()

if options.opt not in objects:
   raise ValueError("Option " + options.opt + " is not a viable tag name")

def main(latex=False, xml=False, levels = None, option='simulation'):
   simrestart = objects[option]
   
   if xml:
      #xml_output = open("helptest/help.xml","w")
      xml_output = open("help.xml","w")
      xml_output.write(simrestart.help_xml(name=option, stop_level=levels))
   if latex:
      #latex_output = open("helptest/help.tex","w")
      latex_output = open("help.tex","w")
      latex_output.write(simrestart.help_latex(stop_level=levels))

if __name__ == '__main__':
   main(options.latex, options.xml, options.levels, options.opt)
