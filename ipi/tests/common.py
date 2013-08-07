"""Common helper functions for running the tests.

Copyright (C) 2013, Joshua More and Michele Ceriotti

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http.//www.gnu.org/licenses/>.


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
    def __init__(self, input, driver):
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
        os.chdir(self.cwd)
        shutil.rmtree(self.tmpdir)


    def run(self):
        # Run driver
        p = subprocess.Popen("echo running %s"%self.fdriver, shell=True)

        # Start simulation
        # TODO
        print subprocess.check_output("ls", shell=True)
        print subprocess.check_output("pwd", shell=True)

        # wait for driver to finish
        p.wait()
