import numpy, math, sys
import cell 
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

def read_array(data):
   """Takes a formatted line with an arra of the form: 
      [array[0], array[1], array[2]], and interprets it"""

   if (data[0] != "[" or data[26] != "," or data[52] != "," or data[78] != "]"):
      print "Error in the array syntax"
   else:
      output = numpy.empty(3)
      length = len(data)
      for i in range(length):
         if data[i] == "D":
            data = data[0:i] + "E" + data[i+1:length]
      try:
         output[0] = float(data[1:26])
         output[1] = float(data[27:52])
         output[2] = float(data[53:78])
         return output
      except ValueError:
         print "Tried to write NaN to array"
         return numpy.zeros(3, float)

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
      return 0.0

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
      
