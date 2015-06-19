import ipi.inputs.simulation as isimulation
from ipi.utils.io.io_xml import xml_parse_file
from ipi.utils.messages import verbosity


def load_simulation(file_name, print_banner=None, print_input=None):

   # TODO: Would this be better as a class method of Simulation?

   # parse the file
   xmlrestart = xml_parse_file(open(file_name))

   # prepare the simulation input object
   input_simulation = isimulation.InputSimulation()

   # check the input and partitions it appropriately
   input_simulation.parse(xmlrestart.fields[0][1])

   # create the simulation object
   simulation = input_simulation.fetch()

   # pipe between the components of the simulation
   simulation.bind()

   # echo the input file if verbose enough
   if verbosity.level > 0:
      print " # i-PI loaded input file: ", file_name
   if verbosity.level > 1:
      print " --- begin input file content ---"
      ifile = open(file_name, "r")
      for line in ifile.readlines():
         print line,
      print " ---  end input file content  ---"
      ifile.close()

   return simulation
