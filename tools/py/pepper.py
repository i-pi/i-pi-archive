#!/usr/bin/env python2

'''
This script cleans python files in-place by calling autopep8 tool.
It searches for python files recursively in the given directory.
Prior to running the script, autopep8 must be installed.

autopep8 is called with the following options:
-i (in-place, modifies files in-place)
-v (verbose)
--select errors (only selected pep8 issues are fixed)
--pep8-passes 2000 (avoids falling into infinite loop by autopep8).

For a full list of issues which can be fixed by autopep8, consult:
https://pypi.python.org/pypi/autopep8

The output of autopep8 is filtered by this script according to verbosity.
There is also option to execute the script only on given files.
In that case, the positional argument is ignored and autopep is run
on the files without checking if they are valid python files.

Syntax:
pepper.py i-pi_root_directory
'''

import argparse
import os
import subprocess
import re
import sys

if __name__ == '__main__':
    # Description of the program after using -h
    parser = argparse.ArgumentParser(
        description='Pepper executes autopep8 '
        'recursive in-place cleaning of python files '
        'according to pep8 style guide in the given directory. '
        'To run it autopep8 must be installed: '
        'https://pypi.python.org/pypi/autopep8  '
        'It recursively looks for all python files in a given directory, '
        'so for example, if you want to clean your i-pi repository, '
        'type: pepper your_ipi_root. '
        'Pepper with recursively search for valid python files and '
        'will clean them. '
        'If you only want to apply it to some files, use --files '
        'option, but BE CAREFUL '
        '- in this mode pepper will not check if they are python files!')
    # There is only one positional argument
    parser.add_argument(
        'path',
        metavar='PATH',
        help='Path to directory to recursively look for python files')
    parser.add_argument('--verbosity',
                        choices=['silent', 'low', 'medium', 'high'],
                        default='medium',
                        help='sets level of verbosity. silent will not print '
                             'anything. low prints messages from this script '
                             'and no autopep8 output, '
                             'medium prints only filenames '
                             'on which script acted and '
                             'high prints everything from autopep8 output')
    parser.add_argument('-f',
                        '--files',
                        type=str,
                        nargs='+',
                        help='run the script only on the given files.'
                             'The positional argument will be ignored. '
                             'pepper will not check if they are python files!')

    args = parser.parse_args()
    path = args.path
    files = args.files
    verbosity = args.verbosity

    # General arguments for pep8
    styles_to_be_corrected = 'E101,E11,E121,E122,E123,E124,E125,E126,E127,E128,' \
                             'E20,E211,E22,E224,E226,E227,E228,E231,E241,E242,E251,E26,E265,E27,' \
                             'E301,E302,E303,E304,E306,' \
                             'E401,' \
                             'E502,' \
                             'W291,W292,W29,' \
                             'W391'
    # Must be written as string since it is passed to command
    number_of_pep8_passes = '2000'
    # in-place, verbose, select only styles which are listed above,
    # pep8-passes limit is to avoid infinite loop
    autopep8_args = ['-i', '-v', '--select', styles_to_be_corrected,
                     '--pep8-passes', number_of_pep8_passes]

    if files is not None:
        # execute autopep8 on these files, without being recursive
        autopep8_args = autopep8_args + files
        if verbosity != 'silent':
            print 'Running autopep8 on the files: ', ' '.join(files)
    else:
        # perform recursive search in the given directory
        if os.path.isdir(path):
            os.chdir(path)
            autopep8_args = autopep8_args + ['-r'] + ['.']
            if verbosity != 'silent':
                print 'Running autopep8 recursively in the directory: ', os.getcwd()
        else:
            print 'The given directory does not exist: ', path
            sys.exit()

    process = subprocess.Popen(
        ['autopep8'] + autopep8_args,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT)
    while process.poll() is None:
        # We must strip, otherwise we get double newline
        line = process.stdout.readline().rstrip()
        if verbosity == 'high':
            print line
        elif verbosity == 'medium':
            # We want to print only filenames that changed
            if re.match('\[file:.*\]', line):
                # Pattern: [file:filename]
                filename_line = re.search(r'\[file:(\S+)\]', line)
                # Print only filename
                print filename_line.group(1)
        # if verbosity is silent or low, do not print output from autopep8

    if verbosity != 'silent':
        print 'autopep8 returned'
