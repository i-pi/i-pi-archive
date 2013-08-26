import os
from setuptools import setup, find_packages 

def read(fname):
   return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
   name='i-PI',
   version='0.9dev',
   description='A Python interface for ab initio path integral molecular dynamics simulations.',
   long_description=read('README.rst'),
   packages=find_packages(),
   author= "Michele Ceriotti",
   author_email = "michele.ceriotti@gmail.com",
   classifiers = [ 'Development Status :: 3 - Alpha' ],
   license='GPLv3'
)
