"""Script which automatically adds the ipi executable to 
/usr/local/pythonX.Y/site-packages.

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
"""

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
