#!/usr/bin/env python2
# pylint: disable=W0601,W0602,W0603,W0604
""" Run regression tests for i-PI. Run in parallel.

Todo:
    * Add support for .gzip files.
    * Improve(Implement) the error handling.

Dependancy:
    * pytest
    * pytest-xdist (required for parallel)
    * numpy

If problem running the test in parallel set processes to 1 in the config file.

The structure of the tree that contains the test must be as follow:

    `-- inputfiles_path    --> Specified in the `config files`
    |-- test1              --> Will be active if present in the `input_files`
    |   |-- something.xml  --> **Only** a single xml file (ipi input)
    |   |-- input_driver   --> Input files for the driver (see below)
    |   `-- output         --> Output files to compare with.
    |-- test_2
    |   |-- something.xml
    |   |-- input_driver
    |   `-- output

The input must contains two specific comment that will tell to specify which
files are needed by the driver and which command should be used to run the
driver. The syntax of the lines is quite strict: only spaces can be different.

To specify which files the driver needs:

```
<!-- The driver will need the following files:  file1, file2, file3 -->
```

To specify the command to call the driver:

```
<!-- Run the driver code with: "i-pi-driver -u -h zundel -m zundel" -->
```

The driver instances will be as many as the number of times the xml contains the
last line.


This script will also try to take care of avoiding socket overlap. The most
most importance think in this regard is that, in case the name of the socket
compare more than once in the "command line" above, only the first one will be
considered. For example: using the command "i-pi-driver -u -m zundel -h zundel"
would not work.

The parallelism is made through pytest. Make sure the plugin xdist is installed!

The path to the config file is hardcoded but can be safly changed.
The words contained in the config file should be listed in `MANDATORY_FIELDS`.


The config file should be commented enough to explain what each keyword does.
Anyway, just to be sure:

[input_files]   => must be followed by the list of test that must be done.
[config]
test_run_path   => path to where the tests will be executed.
inputfiles_path => path to the folder containing all the tests.
processes       => maximum number of concurrent tests that can be executed.
ipi_command     => the command to run i-pi.
initial_address => this is important mostly for inet socket.
precision       => precision (number of decimal) used when comparing numbers.

The output in parallel is very minimal. The advice is to test in serial the
cases that fail in parallel.

To be able to have a kind of debugging output run:
```
pytest -s regtests.py
```

"""

import copy
import os
import re
import shutil
import shlex
import subprocess as sbps
import sys
import time
import xml.etree.ElementTree as etree
from glob import glob

import numpy as np
import numpy.testing as npt
import pytest

# pylint: disable=import-error
from ipi.engine.outputs import CheckpointOutput, PropertyOutput, \
    TrajectoryOutput
from ipi.engine.properties import getkey
from ipi.inputs.simulation import InputSimulation
from ipi.utils.io.inputs import io_xml
# pylint: enable=import-error

######### Hardocoded settings #########
RUNSCRIPT_DIR = os.getcwd()
CONFIG_FILE = './regtests.config'
MANDATORY_FIELDS = [
    'config/inputfiles_path',
    'config/test_run_path',
    'config/ipi_command',
    'config/processes',
    'config/initial_address',
    'config/precision',
    'input_files',
]
#######################################



def main():
    """ Main function: do some check and start the pytest process.

    Pay attention that some operation is done when the module is loaded.
    (see the end of this file)
    """

    if CONFIG['config']['processes'] == 1:
        pytest.main(['-v', os.path.abspath(__file__)])
    else:
        pytest.main(['-v', '-n', CONFIG['config']['processes'],
                     os.path.abspath(__file__)])

def parse_config():
    """ The config file contains the names of the input files and some options.
    """

    # Keyword in the file:
    comment_rgx = re.compile(r'^\s*\#.*')
    section_rgx = re.compile(r'^\s*\[(.*)\]\s*(?:\#.*)?')
    pair_rgx = re.compile(r'^\s*([\w]*)\s*=\s*([\w\'\"\.\"\/\-]*)\s*(?:\#.*)?')
    name_rgx = re.compile(r'\s*\*\s*([\.\w]*)\s*(?:\#.*)?')

    global CONFIG
    global TEST_RUN_PATH
    global INPUTFILES_PATH
    CONFIG = {}

    # Read the config file
    with open(CONFIG_FILE, 'r') as cfg:
        for line in cfg:
            if comment_rgx.match(line):
                continue
            elif section_rgx.match(line):
                section = section_rgx.match(line).group(1)
                # config[section] = None
            elif pair_rgx.match(line):
                key = pair_rgx.match(line).group(1)
                value = pair_rgx.match(line).group(2).strip('"\' ')
                try:
                    CONFIG[section][key] = value
                except (KeyError, AttributeError):
                    CONFIG[section] = {}
                    CONFIG[section][key] = value

            elif name_rgx.match(line):
                try:
                    CONFIG[section].append(name_rgx.match(line).group(1))
                except (KeyError, AttributeError):
                    CONFIG[section] = []
                    CONFIG[section].append(name_rgx.match(line).group(1))

    # Check if all the mandatory config are present
    for fields in [x.split('/') for x in MANDATORY_FIELDS]:
        _tmp_config = copy.deepcopy(CONFIG)
        _old_field = 'root'
        for field in fields:
            try:
                _tmp_config = _tmp_config[field]
                _old_field = field
            except KeyError:
                raise InputError(field, _old_field)


    # Pay attention to avoid the same name/address for different sockets
    filename_list = copy.deepcopy(CONFIG['input_files'])
    initial_address = int(CONFIG['config']['initial_address'])
    for _ii, name in enumerate(filename_list):
        CONFIG['input_files'][_ii] = (name, initial_address)
        initial_address += 1


    TEST_RUN_PATH = os.path.abspath(CONFIG['config']['test_run_path'])
    INPUTFILES_PATH = os.path.abspath(CONFIG['config']['inputfiles_path'])


def general_initialization():
    """ Check that the value given in the config file are usable.

    Do those test that could make something bad happen later...
    """

    # Check if the test_run_path exists and that is a directory
    if os.path.exists(TEST_RUN_PATH) and \
       os.path.isdir(TEST_RUN_PATH):
        answer = 'y'
        while True:
            # answer = raw_input('The directory within you want to run the '+\
            #                    'tests exists.\nDo you want to '+\
            #                    'delete all the contents of\n %s \n[y/n]' \
            #                    % str(TEST_RUN_PATH))
            if answer.lower() == 'n':
                raise KeyboardInterrupt('Usere selected no!')

            elif answer.lower() == 'y':
                shutil.rmtree(TEST_RUN_PATH)
                break
            else:
                sys.stderr.write('The answer must be "y" or "n"')

    sys.stderr.write('Creating the directory %s\n' % TEST_RUN_PATH)
    os.mkdir(TEST_RUN_PATH)


    if len(CONFIG['input_files']) == 0:
        raise RuntimeError('No input files provided')


def initialize_test(test_name): # pylint: disable=too-many-locals
    """ Initialize the tests creating the suitable directories trees.

    All the test will be performed in a directory taken from the config file.
    Each test will have a own subdirectory that will contain other three
    subdirectory.

    `-- TEST_RUN_PATH
    |-- test_path
    |   |-- io           --> All the file produced while running the sofware
    |   |-- input        --> Contains the original data
    |   `-- output       --> Contains only the compared data
    |-- test_2
    |   |-- io
    |   |-- input
    |   `-- output

    This functions prepare the tree folders for a given test.


    Args:
        test_name: A tuple that contains the name of the directory containing
            the test and a number that is used to avoid overlap between sockets.

    Returns:
        The command line string necessary to run the driver code.

    """

    test_name, socket_number = test_name

    driver_command_rgx = re.compile(r'\<\!\-\-\s*[Rr]un\s*the\s*driver\s*code\s*with\:\s*"(.*)"\s*\-\-\>',
                                    re.MULTILINE)
    needed_files_rgx = re.compile(r'\<\!\-\-\s*[Tt]he\s*driver\s*will\s*need\s*the\s*following\s*files\:\s*(.*)\s*\-\-\>',
                                  re.MULTILINE)

    test_path = os.path.join(TEST_RUN_PATH, test_name)
    needed_paths = ['io'] # input is created when copying the input

    # Create the needed directory
    os.mkdir(test_path)
    for path in needed_paths:
        _path = os.path.join(test_path, path)
        os.mkdir(_path)

    # Copy the original files from the INPUTFILES_PATH
    _src = os.path.join(INPUTFILES_PATH, test_name)
    _dst = os.path.join(test_path, 'input')
    shutil.copytree(_src, _dst)

    # Check how many xml are there: if only one assume that is the input,
    #+otherwise take the ipi_input.xml
    xml_file_path_regex = os.path.join(test_path, 'input', '*.xml')
    xml_files = glob(xml_file_path_regex)
    if len(xml_files) == 1:
        xml_file_path = xml_files[0]
    else:
        xml_file_path = os.path.join(test_path, 'input', 'ipi_input.xml')


    # Retrieve the "driver command" and the "needed files" from the xml file
    with open(xml_file_path) as _file:
        _text = _file.read()
    driver_command = [x.group(1) for x in driver_command_rgx.finditer(_text)]
    _needed_files = [x.group(1) for x in needed_files_rgx.finditer(_text)]
    needed_files = [os.path.basename(xml_file_path)]
    for _xx in _needed_files:
        needed_files += _xx.split(',')

    # Copy all the needed files to run the driver and the ipi input into the io
    _dst = os.path.join(test_path, 'io')
    _src_dir = os.path.join(test_path, 'input')
    print 'Files that will be used: ', ', '.join(needed_files)
    print 'Command to run the driver: ', ' | '.join(driver_command)
    for _file in needed_files:
        _src = os.path.join(_src_dir, _file.strip())
        shutil.copy(_src, _dst)


    # Make sure socket will not overlap
    _input_xml = os.path.join(_dst, os.path.basename(xml_file_path))
    xml = etree.parse(_input_xml)
    for ffsocket in xml.findall('./ffsocket'):
        socket_mode = ffsocket.attrib['mode'] # pylint: disable=no-member
        if socket_mode.lower() == 'unix':
            address = ffsocket.find('./address').text # pylint: disable=no-member
            ffsocket.find('./address').text = address+str(socket_number)
            if os.path.exists(os.path.join('/tmp/', 'ipi_'+address+str(socket_number))):
                os.remove(os.path.join('/tmp/', 'ipi_'+address+str(socket_number)))
            for _ii, cmd in enumerate(driver_command):
#                driver_command[_ii].replace(address.strip(), address.strip()+str(socket_number))
                cmd = driver_command[_ii].split()
                for _ww in xrange(len(cmd)):
                    if cmd[_ww].strip() == address.strip():
                        cmd[_ww] = address.strip()+str(socket_number)
                        break
                driver_command[_ii] = ' '.join(cmd)

        else:
            port = ffsocket.find('./port').text
            ffsocket.find('./port').text = str(socket_number)
            for cmd in driver_command:
                cmd.replace(port, str(socket_number))

#    indent(xml)
    xml.write(_input_xml)

    return driver_command, os.path.basename(xml_file_path)


def run_computation(test_name, driver_command, ipi_input_file):
    """ Each test need an instance of i-pi and at least one of driver.

    While the command to run i-PI is quite standard, the command to start the
    driver could be very different for different examples. Thus, the driver
    command is located inside a specific comment in the XML file and retrieved
    during the initialization of the tests.

    """

    # Move to the right directory
    _dir = os.path.join(TEST_RUN_PATH, test_name, 'io')
    ipi_out_path = os.path.join(_dir, 'ipi_output.out')
    driver_out_path = os.path.join(_dir, 'driver_output.out')
    ipi_input_path = os.path.join(_dir, ipi_input_file)
    os.chdir(_dir)

    # Run the i-pi code
    ipi_command = shlex.split(CONFIG['config']['ipi_command'] +\
                              ' ' + ipi_input_path)
    with open(ipi_out_path, 'w') as ipi_out:
        ipi_proc = sbps.Popen(ipi_command,
                              bufsize=0,
                              stdout=ipi_out,
                              stderr=sbps.STDOUT)

    # Sleep few seconds waiting for the ipi server start
    time.sleep(2)

    # Run the driver code
    for cmd in driver_command:
        cmd = shlex.split(cmd)
        with open(driver_out_path, 'w') as driver_out:
            driver_proc = sbps.Popen(cmd,
                                     bufsize=0,
                                     stdin=None,
                                     stdout=driver_out,
                                     stderr=sbps.STDOUT)

    driver_proc.communicate()
    chdir_back()
    return (ipi_proc, driver_proc)


def filesname_to_compare(test_name, input_file):
    """ The test results shold be analyzed number by numbers.

    The idea is that the testing input should never change, then the files
    should be always the same. It would probably be better, anyway, to use the
    i-pi infrastructure to retrieve the right position of the data. In fact,
    this would work as a further testing.

    Returns:
        lprop
        nprop
    """

    # Move to the io directory where all the new files should be find
    test_dir = os.path.join(TEST_RUN_PATH, test_name, 'io')
    orig_dir = os.path.join(TEST_RUN_PATH, test_name, 'input')
    inputfile = os.path.join(test_dir, input_file)
    os.chdir(test_dir)

    # opens & parses the input file
    ifile = open(inputfile, "r")
    xmlrestart = io_xml.xml_parse_file(ifile) # Parses the file.
    ifile.close()

    isimul = InputSimulation()
    isimul.parse(xmlrestart.fields[0][1])

    simul = isimul.fetch()

    # reconstructs the list of the property and trajectory files
    lprop = [] # list of property files
    ltraj = [] # list of trajectory files
    for o in simul.outtemplate:
        if isinstance(o, CheckpointOutput):   # properties and trajectories are output per system
            pass
        elif isinstance(o, PropertyOutput):
            nprop = []
            isys = 0
            for _ in simul.syslist:   # create multiple copies
                filename = o.filename
                nprop.append( { "old_filename" : os.path.join(orig_dir, filename),
                                "new_filename" : os.path.join(test_dir, filename),
                                "stride": o.stride,
                                "properties": o.outlist,} )
                isys+=1
            lprop.append(nprop)

        elif isinstance(o, TrajectoryOutput):   # trajectories are more complex, as some have per-bead output
            if getkey(o.what) in ["positions", "velocities",
                                  "forces", "extras"]:   # multiple beads
                nbeads = simul.syslist[0].beads.nbeads
                for b in range(nbeads):
                    ntraj = []
                    isys=0
                    # zero-padded bead number
                    padb = ( ("%0" + str(int(1 + np.floor(np.log(nbeads)/np.log(10)))) + "d") % (b) )
                    for _ in simul.syslist:
                        if (o.ibead < 0 or o.ibead == b):
                            if getkey(o.what) == "extras":
                                filename = o.filename+"_" + padb
                            else:
                                filename = o.filename+"_" + padb + "." + o.format
                            ntraj.append({ "old_filename" : os.path.join(orig_dir, filename),
                                           "format" : o.format,
                                           "new_filename" : os.path.join(test_dir, filename),
                                           "stride": o.stride,
                                           "what": o.what,})
                        isys += 1
                    if ntraj != []:
                        ltraj.append(ntraj)

            else:
                ntraj=[]
                isys = 0
                for _ in simul.syslist:   # create multiple copies
                    filename=o.filename
                    filename=filename+"."+o.format
                    ntraj.append( { "old_filename" : os.path.join(orig_dir, filename),
                                    "new_filename" : os.path.join(test_dir, filename),
                                    "format" : o.format,
                                    "stride": o.stride,} )

                    isys+=1
                ltraj.append(ntraj)

    chdir_back()
    return ltraj, lprop


def compare_files(test_name, ltraj, lprop):
    test_dir = os.path.join(TEST_RUN_PATH, test_name, 'io')
    os.chdir(test_dir)
    err = False

    for prop in lprop:
        for sprop in prop:
            old_content = np.loadtxt(sprop['old_filename'])
            new_content = np.loadtxt(sprop['new_filename'])

            try:
                npt.assert_array_almost_equal(old_content, new_content,
                                              int(CONFIG['config']['precision']))
            except AssertionError:
                name = os.path.basename(sprop['old_filename'])
                print 'Differences in the %s file' % name
                err = True


    for traj in ltraj:
        for straj in traj:
            new_w_list = []
            old_w_list = []
            name = os.path.basename(straj['old_filename'])
            with open(straj['old_filename']) as old_content:
                with open(straj['new_filename']) as new_content:
                    line_c = 1
                    for old_line, new_line in zip(old_content, new_content):
                        word_c = 1
                        for old_w, new_w in zip(old_line.split(), new_line.split()):
                            # if contains_string.match(old_w):
                            try:
                                old_w_list.append(float(old_w))
                                new_w_list.append(float(new_w))
                            except ValueError:
                                try:
                                    assert old_w == new_w
                                except AssertionError:
                                    print 'Differences at line %d word %d of file  %s' % (line_c, word_c, name)
                                    err = True
                            word_c += 1
                        line_c += 1

                    try:
                        npt.assert_array_almost_equal(np.array(new_w_list), np.array(old_w_list),
                                                      int(CONFIG['config']['precision']))
                    except AssertionError:
                        print 'Differences in the %s file' % name
                        err = True

    if err == True:
        raise AssertionError

    chdir_back()


def chdir_back():
    """ Bring back to the original directory where the script has been called.
    """
    os.chdir(RUNSCRIPT_DIR)





class InputError(Exception):
    """Exception raised for errors in the input.

    Attributes:
        expr -- input expression in which the error occurred
        msg  -- explanation of the error
    """

    def __init__(self, field, oldfield): #pylint: disable=W0231
        self.field = str(field)
        self.old_field = str(oldfield)

    def __str__(self):
        return 'Keyword: %s not found under %s' % \
                                 (self.field, self.old_field)



parse_config()
general_initialization()

@pytest.mark.parametrize("_test", CONFIG['input_files'],
                         ids=[x[0] for x in CONFIG['input_files']])
def test(_test):

    # print '############# Test Name:', _test[0]
    print
    driver_command, ipi_input_file = initialize_test(_test)
    run_computation(_test[0], driver_command, ipi_input_file)

    # Avoid ipi output
    devnull = open('/dev/null', 'w')
    oldstdout_fno = os.dup(sys.stdout.fileno())
    os.dup2(devnull.fileno(), 1)
    lprop, nprop = filesname_to_compare(_test[0], ipi_input_file)
    os.dup2(oldstdout_fno, 1)
    compare_files(_test[0], lprop, nprop)




if __name__ == '__main__':
    main()
