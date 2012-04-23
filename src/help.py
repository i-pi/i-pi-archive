import sys
from engine import simulation
from inputs.simulation import InputSimulation
from utils.io.io_xml import *

def main(latex=False, xml=False):
   simrestart = InputSimulation()
   
   if xml:
      xml_output = open("helptest/help.xml","w")
      xml_output.write(simrestart.help_xml("simulation"))
   if latex:
      latex_output = open("helptest/help.tex","w")
      latex_output.write(simrestart.help_latex())

if __name__ == '__main__':
   main(True, True)
