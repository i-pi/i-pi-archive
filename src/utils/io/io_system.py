import numpy, math, sys, string
import cell, dynamics, rp_dynamics, engine, rp_engine, thermostat, langevin, forces, PILE, barostat, Bussi
import xml.sax.handler, xml.sax, pprint

def print_pdb(atoms, ncell, filedesc = sys.stdout):
   """Takes the system and gives pdb formatted output for the unit cell and the
      atomic positions """

   a, b, c, alpha, beta, gamma = cell.h2abc(ncell.h.get())
   alpha *= 180.0/math.pi #radian to degree conversion
   beta  *= 180.0/math.pi
   gamma *= 180.0/math.pi
   
   z = 1 #number of polymeric chains in a unit cell. I can't decide if 1 or 0 is more sensible for this...

   filedesc.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%s%4i\n" % (a, b, c, alpha, beta, gamma, " P 1        ", z))

   for i in range(0,len(atoms)): 
      filedesc.write("ATOM  %5i %4s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2i\n" % (i+1, atoms[i].name.get(),' ','  1',' ',1,' ',atoms[i].q.get()[0],atoms[i].q.get()[1],atoms[i].q.get()[2],0.0,0.0,'  ',0))

   filedesc.write("END\n")

def print_pdb_RP(systems, filedesc=sys.stdout):

   for syst in systems:
      a, b, c, alpha, beta, gamma = cell.h2abc(syst.cell.h.get())
      alpha *= 180.0/math.pi #radian to degree conversion
      beta  *= 180.0/math.pi
      gamma *= 180.0/math.pi
      
      z = 1 #number of polymeric chains in a unit cell. I can't decide if 1 or 0 is more sensible for this...
   
      filedesc.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%s%4i\n" % (a, b, c, alpha, beta, gamma, " P 1        ", z))
   
      for i in range(0,len(syst.atoms)): 
         filedesc.write("ATOM  %5i %4s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2i\n" % (i+1, syst.atoms[i].name.get(),' ','  1',' ',1,' ',syst.atoms[i].q.get()[0],syst.atoms[i].q.get()[1],syst.atoms[i].q.get()[2],0.0,0.0,'  ',0))

   filedesc.write("END\n")

def read_pdb(filedesc):
   """Takes a pdb-style file and creates a system with the appropriate unit
      cell and atom positions"""

   header = filedesc.readline()
   a = float(header[6:15]);      b = float(header[15:24]);    c = float(header[24:33]);
   alpha = float(header[33:40]); beta = float(header[40:47]); gamma = float(header[47:54]);
   alpha *= math.pi/180.0;       beta *= math.pi/180.0;       gamma *= math.pi/180.0
   cell = numpy.array([a, b, c, alpha, beta, gamma])

   atoms = []
   natoms = 0
   body = filedesc.readline()
   while body != '':
      natoms += 1
      name = body[12:16]
      x = float(body[31:39])
      y = float(body[39:47])
      z = float(body[47:55])
      pos = numpy.array([x, y, z])
      atoms.append([name, pos]) 
      body = filedesc.readline()
   return atoms, cell, natoms

def xml_write(system, namedpipe):
   """Writes an xml-compliant file to file with the atoms positions and cell
      variables"""

   tab = "   "
   namedpipe.write("<?xml version=\"1.0\"?>\n")
   namedpipe.write("<System>\n")
   namedpipe.write(tab + "<natoms>" + str(len(system.atoms)) + "</natoms>\n")

   for i in range(len(system.atoms)):
      atom_q = system.atoms[i].q.get()
      namedpipe.write(tab + "<Atom_vec>\n")
      namedpipe.write(tab + tab + "<q>[" + str(atom_q[0]) + "," + str(atom_q[1]) + "," + str(atom_q[2]) + "]</q>\n")
      namedpipe.write(tab + "</Atom_vec>\n")

   h = system.cell.h.get()
   ih = system.cell.ih.get()
   namedpipe.write(tab + "<Cell_vec>\n")
   namedpipe.write(tab + tab + "<h>[" + str(h[0,0]) + "," + str(h[1,0]) + "," + str(h[2,0]) + "]</h>\n")
   namedpipe.write(tab + "</Cell_vec>\n")
   namedpipe.write(tab + "<Cell_vec>\n" )
   namedpipe.write(tab + tab + "<h>[" + str(h[0,1]) + "," + str(h[1,1]) + "," + str(h[2,1]) + "]</h>\n")
   namedpipe.write(tab + "</Cell_vec>\n")
   namedpipe.write(tab + "<Cell_vec>\n")
   namedpipe.write(tab + tab + "<h>[" + str(h[0,2]) + "," + str(h[1,2]) + "," + str(h[2,2]) + "]</h>\n")
   namedpipe.write(tab + "</Cell_vec>\n")

   namedpipe.write(tab + "<Cell_vec>\n")
   namedpipe.write(tab + tab + "<ih>[" + str(ih[0,0]) + "," + str(ih[1,0]) + "," + str(ih[2,0]) + "]</ih>\n")
   namedpipe.write(tab + "</Cell_vec>\n")
   namedpipe.write(tab + "<Cell_vec>\n" )
   namedpipe.write(tab + tab + "<ih>[" + str(ih[0,1]) + "," + str(ih[1,1]) + "," + str(ih[2,1]) + "]</ih>\n")
   namedpipe.write(tab + "</Cell_vec>\n")
   namedpipe.write(tab + "<Cell_vec>\n")
   namedpipe.write(tab + tab + "<ih>[" + str(ih[0,2]) + "," + str(ih[1,2]) + "," + str(ih[2,2]) + "]</ih>\n")
   namedpipe.write(tab + "</Cell_vec>\n")

   namedpipe.write("</System>\n")

def xml_terminate(namedpipe):
   """Writes a minimal xml-compliant file, which is used to terminate the 
      external program"""

   namedpipe.write("<?xml version=\"1.0\"?>\n")
   namedpipe.write("<terminate></terminate>\n")

class System_read(xml.sax.handler.ContentHandler):
   """Handles reading the xml file containing the force calculations"""

   def __init__(self):
      self.in_pot = False
      self.in_f = False
      self.in_vir = False
      self.pot = ""
      self.f = []
      self.vir = []
      self.buffer = ""

   def startElement(self, name, attributes):
      if name == "pot":
         self.in_pot = True
      elif name == "atom_f":
         self.in_f = True
      elif name == "vir_column":
         self.in_vir = True
      elif name == "f":
         self.buffer = ""
      elif name == "x":
         self.buffer = ""

   def characters(self, data):
      if self.in_pot:
         self.in_pot = False
         self.pot += data
      elif self.in_f:
         self.buffer += data
      elif self.in_vir:
         self.buffer += data

   def endElement(self, name):
      if name == "atom_f":
         self.in_f = False
         self.f.append(self.buffer)
      elif name == "vir_column":
         self.in_vir = False
         self.vir.append(self.buffer)

class Init_object:
   def __init__(self):
      self.kind = ""
      self.parameters = {}

class File_object:
   def __init__(self):
      self.kind = ""
      self.quantity_list = []

class Ensemble_object:
   def __init__(self):
      self.ffield = Init_object()
      self.thermo = Init_object()
      self.cell_thermo = Init_object()
      self.barostat = Init_object()

class Init_read(xml.sax.handler.ContentHandler):
   """Handles reading the xml file containing the initialisation input"""

   def __init__(self):
      self.in_ensemble_kind = False
      self.ensemble_kind = ""
      self.ensemble_object = Ensemble_object()
      self.in_path_integral = False
      self.path_integral = False
      self.in_sys_file = False
      self.sys_file = ""

      self.in_ffield = False
      self.ffield = Init_object()
      self.in_thermo = False
      self.thermo = Init_object()
      self.in_cell_thermo = False
      self.cell_thermo = Init_object()
      self.in_barostat = False
      self.barostat = Init_object()

      self.in_kind = False
      self.in_parameters = False

      self.in_temp = False
      self.temp = 1.0
      self.in_nbeads = False
      self.nbeads = 8

      self.in_Pext = False
      self.Pext = []
      
      self.in_ref_cell = False
      self.ref_cell = []
      
      self.in_x = False

      self.in_time_step = False
      self.time_step = 1.0
      self.in_step_number = False
      self.step_number = 100

      self.in_output_vars = False
      self.output_vars = []
      self.in_output_file = False
      self.output_file = File_object()

      self.in_file_name = False
      self.in_quantity = False

   def startElement(self, name, attributes):
      if name == "ensemble_kind":
         self.ensemble_kind = ""
         self.in_ensemble_kind = True
      elif name == "path_integral":
         self.path_integral = False
         self.in_path_integral = True
      elif name == "sys_file":
         self.sys_file = ""
         self.in_sys_file = True
      elif name == "ffield":
         self.ffield = Init_object()
         self.in_ffield = True
      elif name == "thermo":
         self.thermo = Init_object()
         self.in_thermo = True
      elif name == "cell_thermo":
         self.cell_thermo = Init_object()
         self.in_cell_thermo = True
      elif name == "barostat":
         self.barostat = Init_object()
         self.in_barostat = True
      elif name == "kind":
         self.in_kind = True
      elif name == "parameters":
         self.in_parameters = True
      elif name == "temp":
         self.temp = 1.0
         self.in_temp = True
      elif name == "nbeads":
         self.nbeads = 8
         self.in_nbeads = True
      elif name == "ref_cell": 
         self.in_ref_cell = True
      elif name == "Pext": 
         self.in_Pext = True
      elif name == "x":
         self.in_x = True
      elif name == "time_step":
         self.time_step = 1.0
         self.in_time_step = True
      elif name == "step_number":
         self.step_number = 100
         self.in_step_number = True
      elif name == "output_vars":
         self.output_vars = []
         self.in_output_vars = True
      elif name == "output_file":
         self.output_file = File_object()
         self.in_output_file = True
      elif name == "file_name":
         self.in_file_name = True 
      elif name == "quantity":
         self.in_quantity = True
      elif name == "simulation":
         pass
      else:
         print "Unrecognized opening tag"
         exit()

   def characters(self, data):
      if self.in_ensemble_kind:
         self.ensemble_kind += string.strip(data)
      elif self.in_path_integral:
         self.path_integral = read_bool(data)
      elif self.in_sys_file:
         self.sys_file += string.strip(data)
      elif self.in_kind:
         if self.in_ffield:
            self.ffield.kind += string.strip(data)
         elif self.in_thermo:
            self.thermo.kind += string.strip(data)
         elif self.in_cell_thermo:
            self.cell_thermo.kind += string.strip(data)
         elif self.in_barostat:
            self.barostat.kind += string.strip(data)
         else:
            print "kind tag in the wrong place"
            exit()
      elif self.in_parameters:
         if self.in_ffield:
            self.ffield.parameters = read_dict(data)
         elif self.in_thermo:
            self.thermo.parameters = read_dict(data)
         elif self.in_cell_thermo:
            self.cell_thermo.parameters = read_dict(data)
         elif self.in_barostat:
            self.barostat.parameters = read_dict(data)
         else:
            print "parameter tag in the wrong place"
            exit()
      elif self.in_temp:
         self.temp = read_float(data)
      elif self.in_nbeads:
         self.nbeads = read_int(data)
      elif self.in_x:
         if self.in_ref_cell:
            self.ref_cell.append(read_array(data))
         elif self.in_Pext:
            self.Pext.append(read_array(data))
         else:
            print "x tag in the wrong place"
            exit()
      elif self.in_time_step:
         self.time_step = read_float(data)
      elif self.in_step_number:
         self.step_number = read_int(data)
      elif self.in_file_name:
         if self.in_output_vars and self.in_output_file:
            self.output_file.kind += string.strip(data)
         else:
            print "file_name tag in the wrong place"
            exit()
      elif self.in_quantity:
         if self.in_output_vars and self.in_output_file:
            self.output_file.quantity_list.append(read_dict(data))
         else:
            print "quantity tag in the wrong place"
            exit()

   def endElement(self, name):
      if name == "ensemble_kind":
         self.in_ensemble_kind = False
      elif name == "path_integral":
         self.in_path_integral = False
      elif name == "sys_file":
         self.in_sys_file = False
      elif name == "ffield":
         self.ensemble_object.ffield = self.ffield
         self.in_ffield = False
      elif name == "thermo":
         self.ensemble_object.thermo = self.thermo
         self.in_thermo = False
      elif name == "cell_thermo":
         self.ensemble_object.cell_thermo = self.cell_thermo
         self.in_cell_thermo = False
      elif name == "barostat":
         self.ensemble_object.barostat = self.barostat
         self.in_barostat = False
      elif name == "kind":
         self.in_kind = False
      elif name == "parameters":
         self.in_parameters = False
      elif name == "temp":
         self.in_temp = False
      elif name == "nbeads":
         self.in_nbeads = False
      elif name == "ref_cell": 
         self.in_ref_cell = False
      elif name == "Pext": 
         self.in_Pext = False
      elif name == "x":
         self.in_x = False
      elif name == "time_step":
         self.in_time_step = False
      elif name == "step_number":
         self.in_step_number = False
      elif name == "output_vars":
         self.in_output_vars = False
      elif name == "output_file":
         if self.in_output_vars:
            self.output_vars.append(self.output_file)
         self.in_output_file = False
      elif name == "file_name":
         self.in_file_name = False 
      elif name == "quantity":
         self.in_quantity = False
      elif name == "simulation":
         pass
      else:
         print "Unrecognized closing tag"
         exit()

def read_dict(data):
   """Takes a line with an map of the form:
      {kw[0]: value[0], kw[1]: value[1], ...}, and interprets it"""

   try:
      begin = data.index("{")
      end = data.index("}")
   except:
      print "Error in map syntax"
      exit()
   
   elements = data.count(",") + 1
   if data.count(":") != elements:
      print "Error in map syntax"
      exit()

   length = len(data)
   comma_list = [i for i in range(length) if data[i] == ","]
   colon_list = [i for i in range(length) if data[i] == ":"]

   try:
      if elements == 1:
         output = {}
         kw = data[begin+1:colon_list[0]]
         kw = string.strip(kw)
         value = data[colon_list[0]+1:end]
         output[kw] = string.strip(value)
      else:
         output = {}
         kw = data[begin+1:colon_list[0]]
         kw = string.strip(kw)
         value = data[colon_list[0]+1:comma_list[0]]
         output[kw] = string.strip(value)

         kw = data[comma_list[elements-2]+1:colon_list[elements-1]]
         kw = string.strip(kw)
         value = data[colon_list[elements-1]+1:end]
         output[kw] = string.strip(value)

         for i in range(1, elements-1):
            kw = data[comma_list[i-1]+1:colon_list[i]]
            kw = string.strip(kw)
            value = data[colon_list[i]+1:comma_list[i]]
            output[kw] = string.strip(value)
      return output
   except ValueError:
      print "Tried to read NaN to float in map"
      exit()

def read_array(data):
   """Takes a line with an array of the form: 
      [array[0], array[1], array[2],...], and interprets it"""

   try:
      begin = data.index("[")
      end = data.index("]")
   except ValueError:
      print "Error in array syntax"
      exit()

   elements = data.count(",") + 1
   length = len(data)
   comma_list = [i for i in range(length) if data[i] == ","]
   for i in range(length):
      if data[i] == "D":
         data = data[0:i] + "E" + data[i+1:length]
  
   try:
      if elements == 1:
         output = numpy.zeros(elements)
         output[0] = float(data[begin+1:end])
      else:
         output = numpy.zeros(elements)
         output[0] = float(data[begin+1:comma_list[0]])
         output[elements-1] = float(data[comma_list[elements-2]+1:end])
         for i in range(1,elements-1):
            output[i] = float(data[comma_list[i-1]+1:comma_list[i]])
      return output
   except ValueError:
      print "Tried to write NaN to array"
      exit()

def read_float(data):
   """Takes a formatted line with a double and interprets it"""

   output = 0.0
   length = len(data)
   for i in range(length):
      if data[i] == "D":
         data = data[0:i] + "E" + data[i+1:length]
   try:
      output = float(data)
      return output
   except ValueError:
      print "Tried to write NaN to float"
      exit()

def read_int(data):
   """Takes a formatted line with a double and interprets it"""

   output = 0
   try:
      output = int(data)
      return output
   except ValueError:
      print "Tried to write NaN to int"
      exit()

def read_bool(data):
   """Takes a formatted line with a double and interprets it"""

   output = False
   try:
      output = bool(int(data))
      return output
   except ValueError:
      print "Tried to write NaN to int in bool"
      exit()

def xml_read(namedpipe):
   """Reads an xml-compliant file and gets the potential, forces and virial"""

   parser = xml.sax.make_parser()
   handler = System_read()
   parser.setContentHandler(handler)
   parser.parse(namedpipe)

   pot=read_float(handler.pot)
   f=numpy.zeros(len(handler.f)*3,float)
   vir=numpy.zeros((3,3),float)
   for i in range(len(handler.f)):
      f[3*i:3*(i+1)] = read_array(handler.f[i])
   for i in range(3):
      vir[:,i] = read_array(handler.vir[i])
   return [ pot, f, vir ]
      
def run_from_xml(input_file):
   parser = xml.sax.make_parser()
   handler = Init_read()
   parser.setContentHandler(handler)
   parser.parse(input_file)

   need_ffield = False
   need_thermo = False
   need_cell_thermo = False
   need_barostat = False
   w = 1.0
   Pext = numpy.zeros((3,3))
   h0 = None

   if handler.ensemble_kind == "NVE":
      need_ffield = True
   elif handler.ensemble_kind == "NVT":
      need_ffield = True
      need_thermo = True
   elif handler.ensemble_kind == "NSH" or handler.ensemble_kind == "NPH":
      need_ffield = True
      need_barostat = True
   elif handler.ensemble_kind == "NST" or handler.ensemble_kind == "NPT":
      need_ffield = True
      need_thermo = True
      need_barostat = True
      need_cell_thermo = True
   else:
      print "Unrecognized ensemble"
      exit()
   
   if handler.path_integral:
      effective_temp = handler.temp*handler.nbeads
   else:
      effective_temp = handler.temp

   if need_ffield:
      if handler.ensemble_object.ffield.kind:
         if handler.ensemble_object.ffield.kind == "forcefield":
            ffield = forces.forcefield()
         elif handler.ensemble_object.ffield.kind == "pipeforce":
            pars = handler.ensemble_object.ffield.parameters
            ffield = forces.pipeforce(pars)
         elif handler.ensemble_object.ffield.kind == "rp_pipeforce":
            if handler.path_integral:
               pars = handler.ensemble_object.ffield.parameters
               ffield = forces.rp_pipeforce(pars)
            else:
               print "rp_pipeforce used for classical run"
               exit()
         else:
            print "Unrecognized forcefield"
            exit()
      else:
         print "No forcefield specified, using default"
         ffield = forces.forcefield()
   
   if need_thermo:
      if handler.ensemble_object.thermo.kind:
         if handler.ensemble_object.thermo.kind == "thermostat":
            thermo = thermostat.thermostat(temp = effective_temp, dt = handler.time_step/2.0)
         elif handler.ensemble_object.thermo.kind == "langevin":
            thermo = langevin.langevin(temp = effective_temp, dt = handler.time_step/2.0, tau = read_float(handler.ensemble_object.thermo.parameters["tau"]))
         elif handler.ensemble_object.thermo.kind == "PILE":
            if path_integral:
               thermo = PILE.PILE(temp = effective_temp, dt = handler.time_step/2.0, tau_0 = read_float(handler.ensemble_object.thermo.parameters["tau_0"]))
            else:
               print "PILE thermostat used for classical run"
               exit()
         else:
            "Unrecognized thermostat"
            exit()
      else:
         print "No thermostat specified, using default"
         thermo = thermostat.thermostat(temp = effective_temp, dt = handler.time_step/2.0)
      
   if need_cell_thermo:
      if handler.ensemble_object.cell_thermo.kind:
         if handler.ensemble_object.cell_thermo.kind == "thermostat":
            cell_thermo = thermostat.thermostat(temp = effective_temp, dt = handler.time_step/2.0)
         elif handler.ensemble_object.cell_thermo.kind == "langevin":
            cell_thermo = langevin.langevin(temp = effective_temp, dt = handler.time_step/2.0, tau = read_float(handler.ensemble_object.cell_thermo.parameters["tau"]))
         elif handler.ensemble_object.cell_thermo.kind == "PILE":
            if handler.path_integral:
               cell_thermo = PILE.PILE(temp = effective_temp, dt = handler.time_step/2.0, tau_0 = read_float(handler.ensemble_object.cell_thermo.parameters["tau_0"]))
            else:
               print "PILE thermostat used for classical run"
               exit()
         else:
            "Unrecognized cell thermostat"
            exit()
      else:
         print "No cell thermostat specified, using default"
         cell_thermo = thermostat.thermostat(temp = effective_temp, dt = handler.time_step/2.0)

#TODO barostat temp = effective temp or actual temp?
#TODO make it check at this stage whether ref_cell and Pext are upper triangular
   if need_barostat:
      if handler.Pext != []:
         for i in range(3):
            Pext[:,i] = handler.Pext[i]
      if handler.ref_cell != []:
         h0 = numpy.zeros((3,3))
         for i in range(3):
            h0[:,i] = handler.ref_cell[i]

      if handler.ensemble_object.barostat.kind:
         w = read_float(handler.ensemble_object.barostat.parameters["w"])
         if handler.ensemble_object.barostat.kind == "Bussi_S":
            baro = Bussi.Bussi_S(pext = Pext, dt = handler.time_step/2.0, w = w, temp = effective_temp)
         elif handler.ensemble_object.barostat.kind == "barostat":
            baro = barostat.barostat(pext = Pext, dt = handler.time_step/2.0)
         else:
            print "Unrecognized barostat"
            exit()
      else:
         print "No barostat specified, using default"
         baro = barostat.barostat(pext = Pext, dt = handler.time_step/2.0)
         

#TODO make it possible to initialise without a file, or at least with different file kinds...

   if handler.path_integral:
      f = open(handler.sys_file, "r")
      syst = rp_engine.RP_sys.from_pdbfile(f, nbeads = handler.nbeads, temp = handler.temp, w = w, h0 = h0, pext = Pext)
   else:
      f = open(handler.sys_file, "r")
      syst = engine.System.from_pdbfile(f, w = w, h0 = h0, pext = Pext)

   if handler.ensemble_kind == "NVE":
      if handler.path_integral:
         ensemble = rp_dynamics.rp_nve_ensemble(syst, ffield, dt = handler.time_step)
      else:
         ensemble = dynamics.nve_ensemble(syst, ffield, dt = handler.time_step)
   elif handler.ensemble_kind == "NVT":
      if handler.path_integral:
         ensemble = rp_dynamics.rp_nvt_ensemble(syst, ffield, thermo, temp = handler.temp, dt = handler.time_step)
      else:
         ensemble = dynamics.nvt_ensemble(syst, ffield, thermo, temp = handler.temp, dt = handler.time_step)

#TODO implement the const_P barostat class so that passing barostat actually works.
#TODO implement the rp ensembles

   elif handler.ensemble_kind == "NSH":
      if handler.path_integral:
         ensemble = rp_dynamics.rp_nsh_ensemble(syst, ffield, baro, dt = handler.time_step)
      else:
         ensemble = dynamics.nsh_ensemble(syst, ffield, baro, dt = handler.time_step)
   elif handler.ensemble_kind == "NST":
      if handler.path_integral:
         ensemble = rp_dynamics.rp_nst_ensemble(syst, ffield, baro, thermo, cell_thermo, temp = handler.temp, dt = handler.time_step)
      else:
         ensemble = dynamics.nst_ensemble(syst, ffield, baro, thermo, cell_thermo, temp = handler.temp, dt = handler.time_step)
   elif handler.ensemble_kind == "NPT":
      if handler.path_integral:
         ensemble = rp_dynamics.rp_npt_ensemble(syst, ffield, baro, thermo, cell_thermo, temp = handler.temp, dt = handler.time_step)
      else:
         ensemble = dynamics.npt_ensemble(syst, ffield, baro, thermo, cell_thermo, temp = handler.temp, dt = handler.time_step)




   for output_file in handler.output_vars:
      if output_file.kind == "sys.stdout":
         output_file.kind = sys.stdout
      else:
         output_file.kind = open(output_file.kind, "w")

      for i in range(len(output_file.quantity_list)):
         file_dict = output_file.quantity_list[i]
         for what in file_dict:
            if what == "name":
               name = output_file.quantity_list[i]["name"]
               try:
                  quantity = getattr(ensemble, name)
                  file_dict["name"] = quantity
               except AttributeError:
                  try:
                     quantity = getattr(ensemble.syst, name)
                     file_dict["name"] = quantity
                  except AttributeError:
                     try:
                        quantity = getattr(ensemble.syst.cell, name)
                        file_dict["name"] = quantity
                     except AttributeError:
                        pass
            elif what == "start" or what == "step":
               output_file.quantity_list[i][what] = read_int(output_file.quantity_list[i][what])
            elif what == "average":
               output_file.quantity_list[i][what] = read_bool(output_file.quantity_list[i][what])
            else:
               print "Unknown key word in quantity ", name
               exit()
         try:
            if output_file.quantity_list[i]["average"]:
               output_file.quantity_list[i]["running_total"] = 0.0
         except KeyError:
            pass

#TODO split up this file into separate files for the pdb reading, xml reading, initialisation and put the bit below into a separate function so that the initialisation and the running of the program are in different functions.
   for istep in range(handler.step_number):
      ensemble.step()
      print "Step ", istep + 1, " of ", handler.step_number
      output(ensemble, istep + 1, handler.output_vars)

def output(ensemble, step_count, output_vars):

   for output_file in output_vars:
      filedesc = output_file.kind
      for quantity in output_file.quantity_list:
         if (step_count-quantity["start"])%quantity["step"] == 0 and step_count >= quantity["start"]:
            if quantity["name"] == "system.pdb":
               print_pdb(ensemble.syst.atoms, ensemble.syst.cell, filedesc = filedesc)
            elif quantity["name"] == "RP_system.pdb":
               print_pdb_RP(ensemble.syst.systems, filedesc = filedesc)
            else:
               try:
                  if quantity["average"]:
                     quantity["running_total"] += quantity["name"].get()
                     filedesc.write(quantity["name"].name + ": " + str(quantity["running_total"]/((step_count-quantity["start"])/quantity["step"]+1)) + " ")
                  else:
                     filedesc.write(quantity["name"].name + ": " + str(quantity["name"].get()) + " ")
               except KeyError:
                  filedesc.write(quantity["name"].name + ": " + str(quantity["name"].get()) + " ")

      filedesc.write("\n")
