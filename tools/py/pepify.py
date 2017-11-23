#!/usr/bin/env python2
"""
This script executes autopep8 in the given directory using options:
-r (recursive, searches recursively for all python files)
-i (in-place, modifies files in-place)
-v (verbose)
--select errors (only selected pep8 issues are fixed)
--pep8-passes 2000 (avoids falling into infinite loop by autopep8).
For a full list of issues which can be fixed by autopep8, consult:
https://pypi.python.org/pypi/autopep8

Prior to running the script, autopep8 must be installed.

Syntax:
pepify.py i-pi_root_directory
"""

import argparse
import os
from subprocess import call

if __name__ == '__main__':
    # Description of the program after using -h
    parser = argparse.ArgumentParser(
        description='This script executes autopep8'
        'recursive in-place fixing of python files'
        'according to pep8 in the given directory.'
        'To run it autopep8 must be installed: '
        'https://pypi.python.org/pypi/autopep8')
    # There is only one positional argument
    parser.add_argument(
        'path_to_ipi',
        metavar='PATH_TO_IPI_ROOT',
        help='Path to your i-pi root')
    args = parser.parse_args()
    path = args.path_to_ipi
    if os.path.isdir(path):
        os.chdir(path)
        print 'Running autopep8 in the directory: ', os.getcwd()
        dirs_in_ipi_to_be_pepified = '.'
        styles_to_be_corrected = 'E101,E11,E121,E122,E123,E124,E125,E126,E127,E128' \
                                 'E20,E211,E22,E224,E226,E227,E228,E231,E241,E242,E251,E26,E265,E27' \
                                 'E301,E302,E303,E304,E306' \
                                 'E401' \
                                 'E502' \
                                 'W291,W292,W29' \
                                 'W391'
        # Must be written as string since it is passed to command
        number_of_pep8_passes = '2000'
        call(['autopep8',
              dirs_in_ipi_to_be_pepified,
              '-r',
              '-i',
              '-v',
              '--select',
              styles_to_be_corrected,
              '--pep8-passes',
              number_of_pep8_passes])
    else:
        print 'The given directory does not exist: ', path
