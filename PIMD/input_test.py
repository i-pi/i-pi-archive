import numpy, math, string
from engine import io_system
import xml.sax, pprint

parser = xml.sax.make_parser()
handler = io_system.Init_read()
parser.setContentHandler(handler)
parser.parse("./input_file.xml")

print handler.ensemble_kind
print handler.path_integral
print handler.sys_file
print handler.ensemble_object.ffield.kind
print handler.ensemble_object.ffield.parameters
print handler.ensemble_object.thermo.kind
print handler.ensemble_object.thermo.parameters
print handler.ensemble_object.cell_thermo.kind
print handler.ensemble_object.cell_thermo.parameters
print handler.ensemble_object.barostat.kind
print handler.ensemble_object.barostat.parameters
print handler.temp
print handler.nbeads
print handler.Pext
print handler.ref_cell
print handler.time_step
print handler.step_number
print handler.output_vars[0].kind
print handler.output_vars[0].quantity_list
print handler.output_vars[1].kind
print handler.output_vars[1].quantity_list

io_system.run_from_xml("./input_file.xml")
