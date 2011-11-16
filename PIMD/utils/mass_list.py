import units

class Mass_list:
   def __init__(self):
      self.masses = dict()
      self.masses["   X"] = 1.0
      self.masses["  Si"] = 28.086*units.amu
      self.masses["   H"] = 1.008*units.amu
      self.masses["  H2"] = 2*1.008*units.amu

mlist = Mass_list()

