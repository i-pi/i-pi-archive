from utils.depend import *
from utils.inputvalue import *
from copy import copy
import engine.outputs
import numpy as np
from engine.properties import getkey

__all__=['InputOutputs', 'InputProperties', 'InputTrajectory',
         'InputCheckpoint']

class InputProperties(InputArray):
   """ Simple input class to describe output for properties.

      Storage class for PropertyOutput.
   """

   default_help = """This class deals with the output of one property. """
   default_label = "PROPERTIES"

   attribs=copy(InputArray.attribs)
   attribs["filename"]=(InputValue,{ "dtype" : str, "default": "out"} )
   attribs["stride"]=(InputValue,{ "dtype" : int, "default": 1 } )

   def __init__(self, help=None, dimension=None, units=None, default=None, dtype=None):
      """ Initializes an InputProperties object by just calling the parent
          with appropriate arguments. """

      super(InputProperties,self).__init__(dtype=str, dimension=dimension, default=default, units=units, help=help)

   def fetch(self):
      """ Returns a PropertyOutput object. """

      return engine.outputs.PropertyOutput(self.filename.fetch(), self.stride.fetch(), super(InputProperties,self).fetch())

   def store(self, prop):
      """ Stores a PropertyOutput object. """

      super(InputProperties,self).store(prop.outlist)
      self.stride.store(prop.stride)
      self.filename.store(prop.filename)

class InputTrajectory(InputValue):
   """ Simple input class to describe output for properties.

      Storage class for TrajectoryOutput.
   """

   default_help = """This class defines how one trajectory file should be output. """
   default_label = "TRAJECTORY"

   attribs=copy(InputValue.attribs)
   attribs["filename"]=(InputValue,{ "dtype" : str, "default": "traj"} )
   attribs["stride"]=(InputValue,{ "dtype" : int, "default": 1 } )
   attribs["format"]=(InputValue,{ "dtype" : str, "default": "xyz" } )

   def __init__(self, help=None, dimension=None, units=None, default=None, dtype=None):
      """ Initializes an InputTrajectory object by just calling the parent. """

      super(InputTrajectory,self).__init__(dtype=str, dimension=dimension, default=default, units=units, help=help)

   def fetch(self):
      """ Returns a TrajectoryOutput object. """

      return engine.outputs.TrajectoryOutput(self.filename.fetch(), self.stride.fetch(), super(InputTrajectory,self).fetch(),self.format.fetch())

   def store(self, traj):
      """ Stores a PropertyOutput object. """

      super(InputTrajectory,self).store(traj.what)
      self.stride.store(traj.stride)
      self.filename.store(traj.filename)
      self.format.store(traj.format)

class InputCheckpoint(InputValue):
   """ Simple input class to describe output for properties.

   Storage class for CheckpointOutput.

   Attributes:
      filename: the (base) name for the checkpoint file.
      stride: a checkpoint should be output every <stride> MD steps.
      overwrite: whether checkpoints should be overwritten, or multiple files output.
   """

   default_help = """This class defines how a checkpoint file should be output. """
   default_label = "CHECKPOINT"

   attribs=copy(InputValue.attribs)
   attribs["filename"]=(InputValue,{ "dtype" : str, "default": "restart"} )
   attribs["stride"]=(InputValue,{ "dtype" : int, "default": 1 } )
   attribs["overwrite"]=(InputValue,{ "dtype" : bool, "default": True } )

   def __init__(self, help=None, dimension=None, units=None, default=None, dtype=None):
      """ Initializes an InputTrajectory object by just calling the parent
          with appropriate arguments. """

      super(InputCheckpoint,self).__init__(dtype=int, dimension=dimension, default=default, units=units, help=help)

   def fetch(self):
      """ Returns a CheckpointOutput object. """

      step=super(InputCheckpoint,self).fetch()
      return engine.outputs.CheckpointOutput(self.filename.fetch(), self.stride.fetch(), self.overwrite.fetch(), step=step )

   def parse(self, xml=None, text=""):

      # just a quick hack to allow an empty element
      try:
         super(InputCheckpoint,self).parse(xml,text)
      except:
         self.value=0

   def store(self, chk):
      """ Stores a PropertyOutput object. """

      super(InputCheckpoint,self).store(chk.step)
      self.stride.store(chk.stride)
      self.filename.store(chk.filename)
      self.overwrite.store(chk.overwrite)

class InputOutputs(Input):
   """ List of outputs input class.

   An example of a dynamic input class: a variable number of tags might be present,
   corresponding to different output requests. This allows for instance to print
   multiple property outputs, with different content and/or output frequency.
     """

   attribs = { "prefix" : ( InputValue, { "dtype" : str,
                                          "default"  : "wrap-pi",
                                          "help"     : "A string that will be the pre-pended to each output file name." })
             }

   dynamic = {  "properties" : (InputProperties, { "help" : "Each of the <properties> tags specify how to create a file in which one or more properties are written, one line per frame. " } ),
               "trajectory" : (InputTrajectory, { "help" : "Each of the <trajectory> tags specify how to create a trajectory file, containing a list of per-atom-coordinate properties. " } ),
               "checkpoint" : (InputCheckpoint, { "help" : "Each of the <checkpoint> tags specify how to create a checkpoint file, which can be used to restart a simulation. " } ),
            }

   default_help = """This class defines how properties, trajectories and checkpoints should be output during the simulation.
    May contain zero, one or many instances of <properties>, <trajectory> or <checkpoint> tags, each giving instructions on how
    one output file should be created and managed. """
   default_label = "OUTPUTS"

   def __init__(self, help=None, dimension=None, units = None, default=None):

      # sets default in a "dynamic" way.
      if default is None:
         default = [ engine.outputs.PropertyOutput("wrap-pi.md", 10, [ "time", "step", "conserved", "temperature", "potential", "kinetic_cv" ] ),
                     engine.outputs.TrajectoryOutput("wrap-pi.pos", 100, "positions", "xyz"),
                     engine.outputs.CheckpointOutput("wrap-pi.checkpoint",1000,overwrite=True)          ]
      super(InputOutputs,self).__init__(help, dimension, units, default)
      if not self._default is None:
         self.store(self._default)
         self._explicit = False

   def fetch(self):
      """ Returs a list of the output objects included in this dynamic container. """

      outlist=[ p.fetch() for (n, p) in self.extra ]
      prefix=self.prefix.fetch()
      if not prefix == "":
         for p in outlist: p.filename=prefix+"."+p.filename

      return outlist

   def store(self, plist):
      """ Stores a list of the output objects, creating a sequence of dynamic containers. """

      self.extra=[]

      self.prefix.store("")
      for el in plist:
         if (isinstance(el, engine.outputs.PropertyOutput)):
            ip=InputProperties(); ip.store(el)
            self.extra.append(("properties", ip) )
         if (isinstance(el, engine.outputs.TrajectoryOutput)):
            ip=InputTrajectory(); ip.store(el)
            self.extra.append(("trajectory", ip) )
         if (isinstance(el, engine.outputs.CheckpointOutput)):
            ip=InputCheckpoint(); ip.store(el)
            self.extra.append(("checkpoint", ip) )
