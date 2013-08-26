"""Common helper functions for running the tests.

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


Functions:
   local: Returns local folder of the tests directory.

Classes:
   TestSimulation: Can be used to test that a particular simulation
      will run properly, given an input file and a driver code.
"""

import glob
import os
import subprocess
import shutil
import tempfile

def local(file=None):
    """Returns local folder of the tests directory.

    Args:
        - file: Append file to the local folder
    """
    if file:
        return os.sep.join(__file__.split(os.sep)[:-1]+[file])
    else:
        return os.sep.join(__file__.split(os.sep)[:-1])

class TestSimulation(object):
   """Simple class used to test various aspects of the simulation.

   Can be used to run an example given the location of an xml
   input file and the location of a suitable driver code.

   Attributes:
      finput: The name of the xml input file
      folder_input: A string giving the directory the input file is held in.
      fdriver: The location of a driver code.
      cwd: Current working directory.
      tmpdir: A temporary directory to run the simulation in.
   """

   def __init__(self, input, driver):
      """Initializes TestSimulation.

      Args:
         input: The name of the xml input file.
         driver: The location of the driver code.
      """

      self.finput = input
      self.folder_input = os.sep.join(input.split(os.sep)[:-1])
      self.fdriver = driver
      self.cwd = os.getcwd()
      self.tmpdir = tempfile.mkdtemp()

      # Copy needed files to tmpdir
      for src in glob.glob("%s/*"%self.folder_input):
          shutil.copy(src, self.tmpdir)

      os.chdir(self.tmpdir)

   def __del__(self):
      """Cleans the temporary directory once the simulation is over."""

      os.chdir(self.cwd)
      shutil.rmtree(self.tmpdir)

   def run(self):
      """Runs the simulation."""

      # Run driver
      p = subprocess.Popen("echo running %s"%self.fdriver, shell=True)

      # Start simulation
      # TODO
      print subprocess.check_output("ls", shell=True)
      print subprocess.check_output("pwd", shell=True)

      # wait for driver to finish
      p.wait()
