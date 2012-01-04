import numpy as np
from utils.io.io_pdb import print_pdb
from utils.io.io_xml import *
from utils.restart import *
from utils.depend import *
import engine.simulation
import sys

class RestartOutput(Restart):
   fields={"stdout_stride": (RestartValue, (int, 100)), "energy_stride": (RestartValue, (int, 1)), 
           "traj_stride": (RestartValue, (int, 500)), "checkpoint_stride": (RestartValue, (int, 1000)),
           "file_prefix": (RestartValue, (str, "")), "quantity_list": (RestartValue, (str, ""))}

   def store(self, output):
      self.stdout_stride.store(output.stdout)
      self.energy_stride.store(output.energy)
      self.traj_stride.store(output.traj)
      self.checkpoint_stride.store(output.check)
      self.file_prefix.store(output.prefix)
      self.quantity_list.store(output.quantity_list)

   def fetch(self):
      return Output(self.stdout_stride.fetch(), self.energy_stride.fetch(), self.traj_stride.fetch(), self.checkpoint_stride.fetch(), self.file_prefix.fetch(), self.quantity_list.fetch())

class Output(dobject):
   def __init__(self, stdout = 100, energy = 1, traj = 500, checkpoint = 1000, prefix = "test", quantity_list = ""):
      self.stdout = stdout
      self.energy = energy
      self.traj = traj
      self.check = checkpoint
      self.prefix = prefix

      self.quantity_list = quantity_list
      self.quantities = self.quantity_list.split(" ")

      self.energy_file = open(self.prefix + ".out", "a")
      self.traj_file = open(self.prefix + ".pdb", "a")

      self.number = 0
      self.file_check()
      self.restart = engine.simulation.RestartSimulation()

   def file_check(self):
      new = False
      self.number += 1
      while (not new):
         try:
            self.check_file = open(self.prefix + ".restart" + str(self.number), "r")
            self.check_file.close()
            self.number += 1
         except IOError:
            self.check_file = open(self.prefix + ".restart" + str(self.number), "w")
            new = True

   def bind(self, properties, simul):
      self.properties = properties
      self.simul = simul
            
   def step(self, step_no, time):
      if not hasattr(self, "write_list"):
         self.write_list = []
         try:
            quant_str="%16s"%("time")
            for quantity in self.quantities:
               if quantity != "":
                  if quantity == "cell_parameters":
                     quant_str += "%16s"%("a") + " " + "%16s"%("b") + " " + "%16s"%("c") + " " + "%16s"%("alpha") + " " + "%16s"%("beta") + " " + "%16s"%("gamma") + " "
                  else:
                     quant_str += "%16s"%(quantity) + " "
                  if quantity != "all":
                     self.write_list.append(self.properties.property_dict[quantity])
                  else:
                     quant_str = "%16s"%("time")
                     for prop in self.properties.property_dict:
                        if prop == "cell_parameters":
                           quant_str += "%16s"%("a") + " " + "%16s"%("b") + " " + "%16s"%("c") + " " + "%16s"%("alpha") + " " + "%16s"%("beta") + " " + "%16s"%("gamma") + " "
                        else:
                           quant_str += "%16s"%(prop) + " "
                        self.write_list.append(self.properties.property_dict[prop])
                     break
            quant_str += "\n"
            self.energy_file.write("# " + quant_str)
         except KeyError:
            print "Quantity ", quantity, " is not available for calculation"
            print "\nAvailable quantities: \n"
            for prop in self.properties.property_dict:
               print prop + " "
            print "\n"
            exit()

      if step_no%self.energy == 0:
         self.energy_file.write("  ")
         self.energy_file.write(write_type(float,time))
         for i in range(len(self.write_list)):
            quantity = self.write_list[i].get()
            try:
               for i in range(len(quantity)):
                  pass
               quantity = list(quantity)
               for i in range(len(quantity)):
                  self.energy_file.write(write_type(float, quantity[i]) + " ")
            except TypeError:
               quantity = float(quantity)
               quantity = write_type(float, quantity)
               self.energy_file.write(quantity + " ")
         self.energy_file.write("\n")

      if step_no%self.traj == 0:
         print_pdb(self.simul.atoms, self.simul.cell, self.traj_file)

      if step_no%self.check == 0 and step_no != 0:
         self.restart.store(self.simul)
         self.check_file.write(self.restart.write("simulation"))
         self.check_file.close()
         self.file_check()

      if step_no%self.stdout == 0:
#TODO put some more relevant stuff here
         sys.stdout.write("step = " + str(step_no) + ", econs = " + str(self.properties.property_dict["econs"].get()) + "\n")
