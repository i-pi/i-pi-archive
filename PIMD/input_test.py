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
print handler.ffield.kind
print handler.ffield.parameters
print handler.thermo.kind
print handler.thermo.parameters
print handler.cell_thermo.kind
print handler.cell_thermo.parameters
print handler.temp
print handler.nbeads
print handler.Pext
print handler.ref_cell
print handler.time_step
print handler.thermostating_steps
print handler.step_number
print handler.output_vars[0].kind
print handler.output_vars[0].parameters
print handler.output_vars[1].kind
print handler.output_vars[1].parameters
