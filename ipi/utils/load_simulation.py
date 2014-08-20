import ipi.inputs.simulation as isimulation
from ipi.utils.io.io_xml import xml_parse_file


def load_simulation(file_name, print_banner=None, print_input=None):

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

   return simulation
