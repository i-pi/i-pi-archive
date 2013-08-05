"""Common helper functions for running the tests.
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
