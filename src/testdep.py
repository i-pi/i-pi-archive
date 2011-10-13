from utils.depend import *

class test2(object):
   # two properties which do not need to be computed
   @depend()
   def velocity(self): pass

   @depend()
   def position(self): pass
   
   #dependent properties
   @depend(0.0,[position])
   def potential(self):
      print "    getting potential"
      return self.position**4
   
   @depend(0.0,[velocity])
   def kinetic(self):
      print "    getting kinetic"
      return self.velocity**2
      
   @depend(3,[potential,kinetic])
   def energy(self):
      print "    getting energy";
      return self.kinetic+self.potential

   @depend(1.0)
   def half(self):
      print "     getting half"
      return self.twice/2

   @depend(2.0,[],[half])
   def twice(self):
      print "     getting twice"
      return self.half*2
   
dp=test2()
dp.position=0.1
dp.velocity=1.0
print "Getting total energy (first evaluation): "
print dp.energy
print "Getting total energy (second evaluation): "
print dp.energy

print "Hard-coded setting of kinetic energy to 4"
dp.kinetic=4.0
print "Getting total energy: "
print dp.energy

print "Now we change the position to 3"
dp.position=3.0
print "Getting total energy: "
print dp.energy

print "Now testing synced objects!"
print "   half: \n", dp.half, "\n   twice: \n",dp.twice

print "Setting half to 10"
dp.half=10.0
print "   half: \n", dp.half, "\n   twice: \n",dp.twice

print "Setting twice to 5"
dp.twice=5.0
print "   half: \n", dp.half, "\n   twice: \n",dp.twice

print "Setting both twice and half to 1"
dp.half=1.0
dp.twice=1.0
print "Result depends on whether we first set half"
print "   half: \n", dp.half, "\n   twice: \n",dp.twice

dp.twice=1.0
dp.half=1.0
print "...or twice"
print "   half: \n", dp.half, "\n   twice: \n",dp.twice

