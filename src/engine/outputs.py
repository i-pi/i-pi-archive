""" Mini classes to deal with output of simulation data. """

import numpy as np
from utils.depend import *
from utils.io.io_xml import *

__all__=[ 'PropertyOutput', 'TrajectoryOutput', 'CheckpointOutput' ]

class PropertyOutput(dobject):
   """ Mini-class dealing with outputting a set of properties to file.

      Does not do any calculation, just manages opening a file, getting data
      from a Properties object and outputting with the desired stride.
   """

   def __init__(self, filename="out", stride=1, outlist=np.zeros(0,np.dtype('|S12')) ):
      """ Initializes a property output stream opening the corresponding file name.
         Also writes out headers.  """
      self.filename=filename
      self.outlist=np.asarray(outlist,np.dtype('|S12'))
      self.stride=stride
      self.out=None


   def bind(self, simul):
      """ Binds output proxy to simulation object. """

      self.simul=simul

      # Checks as soon as possible if some asked-for properties are missing or mispelled
      for what in self.outlist:
         if '(' in what:
            argstart = what.find('(')
            key = what[0:argstart]
            if not key in self.simul.properties.property_dict.keys():
               print "Computable properties list: ", self.properties.property_dict.keys()
               raise KeyError(key + " is not a recognized property")

         elif not what in self.simul.properties.property_dict.keys():
            print "Computable properties list: ", self.properties.property_dict.keys()
            raise KeyError(what + " is not a recognized property")

      self.open_stream()

   def open_stream(self):
      """ Opens the output stream. """

      try:
         self.out=open(self.filename, "a")
      except:
         raise ValueError("Could not open file "+self.filename+" for output")

      ohead = "# "
      for l in self.outlist:
         ohead += "%16s"%(l) + " "
      self.out.write(ohead + "\n")


   def close_stream():
      """ Closes the output stream. """

      self.out.close()

   def write(self):
      """Outputs the required properties of the system.

      Note that properties are outputted using the same format as for the
      output to the xml checkpoint files, as specified in io_xml.

      Raises:
         KeyError: Raised if one of the properties specified in the output list
            are not contained in the property_dict member of properties.
      """


      if not (self.simul.step+1) % self.stride == 0:   return
      self.out.write("  ")
      for what in self.outlist:
         try:
            quantity = self.simul.properties[what]
         except KeyError:
            raise KeyError(what + " is not a recognized property")
         if not hasattr(quantity,"__len__") :
            self.out.write(write_type(float, quantity) + "   ")
         else:
            self.out.write(" [ ")
            for el in quantity:
               self.fout.write(write_type(float, el) + " ")
            self.out.write(" ] ")

      self.out.write("\n")
      self.out.flush()


class TrajectoryOutput(dobject):
   """ Mini-class dealing with outputting atom-based properties as a trajectory file.

      Does not do any calculation, just manages opening a file, getting data
      from a Trajectories object and outputting with the desired stride.
   """


   def __init__(self, filename="out", stride=1, what="", format="xyz"):
      """ Initializes a property output stream opening the corresponding file name.
         Also writes out headers.  """
      self.filename=filename
      self.what=what
      self.stride=stride
      self.format=format
      self.out=None


   def bind(self, simul):
      """ Binds output proxy to simulation object. """

      self.simul=simul
      self.open_stream()

   def open_stream(self):
      """ Opens the output stream(s). """

      if self.what in [ "positions", "velocities", "forces" ]:
         # must write out trajectories for each bead, so must create b streams
         self.out=[]
         for b in range(self.simul.beads.nbeads):
            padb=( ("%0"+str(int(1+np.floor(np.log(self.simul.beads.nbeads)/np.log(10))))+"d") % (b) )
            try:
               self.out.append( open(self.filename+"_"+padb+"."+self.format, "a") )
            except:
               raise ValueError("Could not open file "+self.filename+"_"+padb+"."+self.format+" for output")

   def close_stream():
      """ Closes the output stream. """

      if hasattr(self.out, "__getitem__"):
         for o in self.out: o.close()
      else:
         self.out.close()

   def write(self):
      """Writes out the required trajectories."""

      if not (self.simul.step+1) % self.stride == 0:   return

      # quick-and-dirty way to check if a trajectory is "global" or per-bead
      # Checks to see if there is a list of files or just a single file.
      if hasattr(self.out, "__getitem__"):
         for b in range(len(self.out)):
            self.simul.trajs.print_traj(self.what, self.out[b], b, format=self.format)
      else:
         self.simul.trajs.print_traj(self.what, self.out, b=0, format=self.format)


class CheckpointOutput(dobject):
   """ Mini-class dealing with outputting checkpoints.

   Saves the complete status of the simulation at regular intervals.

   Attributes:
      filename: the (base) filename for the checkpoint file file.
      step: the number of times a checkpoint has been written out.
      stride: will print out a checkpoint every <stride> MD steps.
      overwrite: if True, the checkpoint file is overwritten at each output.
         if False, will output to filename_step. note that no check is done
         on whether filename_step exists already.
   """


   def __init__(self, filename="restart", stride=1000, overwrite=True, step=0):
      """ Initializes a checkpoint output proxy. """

      self.filename=filename
      self.step=step
      self.stride=stride
      self.overwrite=overwrite


   def bind(self, simul):
      """ Binds output proxy to simulation object. """

      import inputs.simulation
      self.simul=simul
      self.status=inputs.simulation.InputSimulation(); self.status.store(simul)

   def store(self):
      self.status.store(self.simul)

   def write(self):
      """Writes out the required trajectories."""

      if not (self.simul.step+1) % self.stride == 0:   return

      if self.overwrite: filename=self.filename
      else: filename=self.filename+"_"+str(self.step)

      self.step+=1    # advances the step counter before saving, so next time the correct index will be loaded.
      self.store()
      check_file = open(filename, "w")
      check_file.write(self.status.write(name="simulation"))
      check_file.close()



