"""Help script which automatically generates help files.

Copyright (C) 2013, Joshua More and Michele Ceriotti

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


This takes an output class specified by the user, and then uses the
automatic help generation functions to generate appropriate help files for this
class, giving information about the tags and the attributes of this class.

A full help message can be found by running 'python help.py -h' or
'python help.py --help'.

Functions:
   help: Writes the help file.
"""

import sys

src_dir = ".."

sys.path.append(src_dir)

from ipi.engine.properties import *
from ipi.utils.io.io_xml import *
from optparse import OptionParser

__all__ = ['help_list', 'list_objects']

list_objects = { 'property_list': Properties(), 
            'trajectory_list': Trajectories()}

usage = "usage: python %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option("-o", action="store", dest="prefix", help="Prefix for the output files", default="help")
parser.add_option("-i", action="store", dest="opt", help="Root object for the help files. Options: ['property_list', 'trajectory_list']", default='property_list')
parser.add_option("-r", action="store_true", dest = "ref", default=False, help="If false, this creates a stand-alone document.")
(options, args) = parser.parse_args()

if options.opt not in list_objects:
   raise ValueError("Option " + options.opt + " is not a viable tag name")

def help_list(option='property_list', prefix="help", standalone=True):
   """Writes the help file.

   Will write a latex file 'prefix.tex'. Can also print out 
   sections of latex documents rather than entire 
   documents, so that we can input them into other latex documents, such as
   the manual.

   Args:
      option: A string specifying which object will be used as the root object
         for the latex and xml files. Defaults to 'property_list'.
      prefix: File prefix for the output files. Defaults to 'help'.
      standalone: Boolean specifying whether the latex file will be a stand-alone
         document, or will instead be intended to be used in a larger document
         with cross references between the different objects.
   """

   simrestart = list_objects[option]
   if option == "property_list":
      idict = list_objects[option].property_dict
   elif option == "trajectory_list":
      idict = list_objects[option].traj_dict
   else:
      raise ValueError("Incorrect option specified.")
   
   latex_output = open(prefix + ".tex","w")
   latex_output.write(help_latex(idict, standalone=standalone))

if __name__ == '__main__':
   help_list(options.opt, options.prefix, not options.ref)
