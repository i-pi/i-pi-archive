import numpy
import math
import thermostat

class Thermo_test(thermostat.Thermostat):     
   
   def __init__(self):
      thermostat.Thermostat.__init__(self)
      self.__tau=1.0
      self.__T=1.0
      self.__S=1.0

