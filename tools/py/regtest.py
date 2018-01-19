#!/usr/bin/env python2

import sys
import re
import time
import os
import shutil
from collections import deque
import xml.etree.ElementTree as etree
import shlex
import subprocess
from ipi.inputs.simulation import InputSimulation
from ipi.utils.io.inputs import io_xml
from ipi.engine.outputs import CheckpointOutput, PropertyOutput, TrajectoryOutput
from ipi.engine.properties import getkey
import ipi.utils.softexit
import numpy as np
import numpy.testing as npt
import argparse

def main():
    arguments = _parser()
    run_dir = os.path.abspath(arguments['run_directory'])
    test_dir = os.path.abspath(arguments['test_cases_directory'])
    is_in_reference_mode = arguments['is_in_reference_mode']
    #force_remove = arguments['force_remove']
    try:
        os.makedirs(run_dir)
    except OSError:
        print "The directory %s exists. Do you want to delete its contents and continue? (y/n)" % run_dir
        if answer_is_y():
            shutil.rmtree(run_dir)
            os.makedirs(run_dir)
        else:
            quit("Terminated")
    if is_in_reference_mode:
        print "Do you agree to replace references if they exist? (y/n)"
        if answer_is_y():
            replace_references = True
        else:
            replace_references = False
        
    list_of_test_candidates = build_test_index(test_dir)
    
    list_of_test_cases = []
    for test_candidate in list_of_test_candidates:
        try:
            list_of_test_cases.append(TestCase(test_candidate))
        except TypeError, e:
            print "Could not create test case:", test_candidate.name
            print str(e)
            continue
    queue_of_test_instances = deque([])
    counter = Counter()
    print
    print "List of test cases to be executed:"
    for test_case in list_of_test_cases:
        test_instance_dir = os.path.join(run_dir,test_case.name)
        queue_of_test_instances.append(TestInstance(test_case,test_instance_dir,counter))
        print  '-->', test_case.name
    queue_of_output = deque([])
    print
    print "Executing following test instances:"
    while queue_of_test_instances:
        test_instance = queue_of_test_instances.popleft()
        try:
            #print '-->', test_instance.name
            sys.stdout.write('--> ' + test_instance.name)
            sys.stdout.flush()
            #print 'test output directory:', test_instance.run_dir
            test_instance.run()
            if not is_in_reference_mode:
                try:
                    test_passed = True
                    differences = test_instance.compare_output_with_reference()
                    for results in differences:
                        if not results.files_are_equal():
                            test_passed = False
                    if test_passed:
                        print '    PASSED'
                    else:
                        print '    FAILED'
                        results.print_differences()
                except ValueError as e:
                    print str(e)
            else:
                try:
                    test_instance.put_output_as_reference()
                except ValueError:
                    if replace_references:
                        test_instance.test_case.remove_reference()
                        test_instance.put_output_as_reference()
                        print '    SUCCESS: References replaced'
                    else:
                        print '    SKIPPING: References not copied'
                
        except (RegtestTimeout, WrongDriverCommand, IPIError, OutputCorrupted) as e:
            print '    ERROR'
            print str(e), 'in:', test_instance.run_dir

            continue
        except WrongIPICommand as e:
            sys.exit(str(e))


def _parser():
    """ Parse the argument lists given as input.

    Return:
        A dictionary containing all the input options and argument.
    """
    parser = argparse.ArgumentParser(description='Recursively looks for '
    'i-pi regression tests '
    'in TEST_CASES_DIRECTORY (default: the current directory), '
    'runs them in RUN_DIRECTORY (default: ./regtest-run) '
    'and compares the output of each test case with the reference output. '
    'With the --create-reference flag the reference outputs are created.')
#    parser.add_argument('--driver-maxtime',
#                        action='store',
#                        type=int,
#                        help=('Wall time for the driver: this is useful when'
#                              'the driver is stuck. The default value is chosen'
#                              'based on the provided regtests.'),
#                        default=Parameters.default_driver_timeout,
#                        dest='driver_maxtime')
#    parser.add_argument('--test',
#                        action='append',
#                        type=str,
#                        help=('Mark a test to be ran. A single test can be '
#                              'specified with "group:test". This option can be '
#                              'used as many times as you want.'),
#                        default=[],
#                        dest='test_list')
    parser.add_argument('--test-cases-directory',
                        action='store',
                        type=str,
                        default=Parameters.default_test_cases_directory,
                        help=('The directory which is recursively searched for regression tests'),
                        dest='test_cases_directory')
    parser.add_argument('--run-directory',
                        action='store',
                        type=str,
                        default=Parameters.default_run_directory,
                        help=('Root directory where the test will be run.'),
                        dest='run_directory')
#    parser.add_argument('--number-of-first-port',
#                        action='store',
#                        type=int,
#                        default=Parameters.default_driver_timeout,
#                        help=('Number of the first port in case of inet '
#                              'socket.'),
#                        dest='number_of_first_port')
    parser.add_argument('--create-reference',
                        action='store_true',
                        default=False,
                        help=('Do not compare output and create the relevant '
                              'reference outputs for test cases.'),
                        dest='is_in_reference_mode')
#    parser.add_argument('--precision',
#                        action='store',
#                        default=Parameters.default_precision,
#                        type=int,
#                        help=('Define the number of decimal used in comparing'
#                              'float.'),
#                        dest='precision')
#    parser.add_argument('-f',
#                        '--force-remove',
#                        action='store_true',
#                        default=False,
#                        help=('Remove run directory and, if in create-reference mode, references without asking the user'),
#                        dest='force_remove')

    return vars(parser.parse_args())

def file_is_test(path_of_xml_file):
    """ Check if an xml file can be used as regtest.
    """

    with open(path_of_xml_file) as _file:
        _text = _file.read()
    regtest_string_rgx = Parameters.compile_regtest_string()
    return (len([x.group(1) for x in regtest_string_rgx.finditer(_text)]) > 0)

def parse_regtest_string(xml_input_path):
    """ Retrieve the commands and the dependecies from the xml input.

    This function parses the comment at the beginning of xml file
    and returns the commands to be used as driver and
    a list of dependencies (files required).

    Args:
        xml_input_path: Path of the xml file.

    Returns:
        commands: list of command to be used to run the driver.
        dependencies: list of files needed to run the test
    """

    regtest_string_rgx = Parameters.compile_regtest_string()
    command_string_rgx = Parameters.compile_command_string()
    dependency_string_rgx = Parameters.compile_dependencies_string()

    with open(xml_input_path) as _buffer:
        _text = _buffer.read()

    regtest_string = []
    dependencies = []
    commands = []

    regtest_string = regtest_string_rgx.findall(_text)

    for _xx in dependency_string_rgx.finditer('\n'.join(regtest_string)):
        dependencies += _xx.group(1).split()
    for _xx in command_string_rgx.finditer('\n'.join(regtest_string)):
        try:
            commands += [_xx.group(2)] * int(_xx.group(1))
        except ValueError:
            commands += [_xx.group(2)]

    return commands, dependencies

def build_test_index(root_path):
    """ Look for a valid xml regtests recursively in the root_path.

    If a directory contains more than a single xml,
    it is skipped and a warning is printed.
    If a directory does not contain any xml, it is skipped without message.

    Args:
        root_path: Parent folder of all the tests.

    Returns:
        list of TestCandidates

    """

    abs_root_path = os.path.abspath(root_path)
    test_list = []

    if not os.path.exists(abs_root_path) or not os.path.isdir(abs_root_path):
        raise RuntimeError('Folder %s does not exist!' % abs_root_path)

    for root, dirs, files in os.walk(abs_root_path):
        test_name = os.path.relpath(root, root_path)
        if test_name == "." : # if you require as test-case path a path that contains a valid .xml file, uses this special name
            test_name = "_root_test_"
        xml_files = [x for x in files if x.endswith('.xml')]
        if len(xml_files) > 1:
            #print ""
            #sys.stderr.write('!W! Skipping test %s: too many xml files!\n'
            #                 % test_name)
            #print "Skipping directory, too many xml files:", root
            continue
        elif len(xml_files) == 1:
            try:
                test_list.append(TestCandidate(test_name,
                              os.path.join(root, xml_files[0])))
            except ValueError:
                continue
    return test_list

def prepare_dir(dir_path):
    """ Creates directory or deletes content if it is not empty.

    Args:
        dir_path: The path of the directory that should be prepared.

    """
    if os.path.exists(dir_path):
        try:
            if os.path.isdir(dir_path):
                shutil.rmtree(dir_path)
            elif os.path.isfile(dir_path):
                os.remove(dir_path)
            else:
                raise RuntimeError
        except:
            raise RuntimeError('I cannot remove the file. '
                               'Try manually and restart this script!')
    os.makedirs(dir_path)
    return

def check_presence_of_dependencies(xml_file, path):
    commands, dependencies = parse_regtest_string(xml_file)
    is_present = True
    #print dependencies
    abs_dependencies = [os.path.join(os.path.abspath(path),x) for x in dependencies]
    #print abs_dependencies
    for files in dependencies:
        file_path = os.path.join(path, files)
        if not os.path.isfile(file_path):
            #print "Missing file", file_path
            is_present = False
    return is_present
    
class TestCase:
    def __init__(self, test_candidate):
        dependencies_dir = os.path.dirname(test_candidate.input_path)
        commands, dependencies = parse_regtest_string(test_candidate.input_path)
        dependencies_list = [os.path.join(dependencies_dir, x) for x in dependencies]
        if not all ([os.path.exists(x) for x in dependencies_list]):
            raise TypeError('Dependencies for file: ' + str(test_candidate.input_path) +
            ' absent in: ' + str(dependencies_dir))
        self.input_path = test_candidate.input_path
        self.name = test_candidate.name
        self.root = os.path.dirname(test_candidate.input_path)
        self.dependencies_list = dependencies_list
        self.dependencies_dir = dependencies_dir
        self.reference_dir = os.path.join(self.root, Parameters.reference_directory)
        
    def how_many_sockets(self):
        xml = etree.parse(self.input_path)
        sockets = 0
        for ffsocket in xml.findall('./ffsocket'):
            sockets += 1
        return sockets
    
    def get_reference_output(self):
        return TestOutput(self.input_path, self.reference_dir)
        
    def remove_reference(self):
        shutil.rmtree(self.reference_dir)

class Counter:
    socket = 10001
    def attribute_socket(self):
        socket = Counter.socket
        Counter.socket += 1
        return socket

class TestInstance:
    def __init__(self, test_case, dir, counter):
        run_dir = os.path.abspath(dir)
        try:
            os.makedirs(run_dir)
        except OSError:
            if os.path.exists(run_dir):
                raise ValueError ("Directory %s exists. The dir parameter must be a path to nonexistent directory" % str(run_dir))
            else:
                raise
        test_dependencies = []
        #print "run_dir", run_dir
        for files in test_case.dependencies_list:
            shutil.copy(files, run_dir)
            test_dependencies.append(os.path.join(run_dir, os.path.basename(files)))
        #print test_dependencies
            
        xml = etree.parse(test_case.input_path)
        commands, dependencies = parse_regtest_string(test_case.input_path)
        for ffsocket in xml.findall('./ffsocket'):
            address_index = counter.attribute_socket()
            socket_mode = ffsocket.attrib['mode']
            if socket_mode.lower() == 'unix':
                address = ffsocket.find('./address').text.strip()
                new_address = ' %s%i ' % (address, address_index)
                ffsocket.find('./address').text = new_address
            else:
                address = ffsocket.find('./port').text.strip()
                new_address = '%i' % address_index
                ffsocket.find('./port').text = new_address

            # Change the address used by the driver too!
            for _ii, _buffer in enumerate(commands):
                # Determine if the address is in a file
                for _word in _buffer.split():
                    _word = os.path.join(run_dir, _word)
                    if os.path.exists(_word) and os.path.isfile(_word):
                        inplace_change(_word, address, new_address)
                        #print _word, address, new_address

                # Replace only the first occurrence of 'address' in the
                #+driver command!
                commands[_ii] = _buffer.replace(address, new_address, 1)


        xml.write(os.path.join(run_dir, 'input.xml'))
        driver_signatures = [] #tuple: command + dir
        for cmd in commands:
            instance_folder = os.path.join(run_dir,
                                           'driver-%i' % len(driver_signatures))
            os.makedirs(instance_folder)
            for _file in test_dependencies:
                shutil.copy(_file, instance_folder)
                #print _file, instance_folder

            driver_command = shlex.split(cmd)
            driver_signatures.append((driver_command, instance_folder))
        self.test_case = test_case
        self.path_ipi_input = os.path.join(run_dir, 'input.xml')
        self.name = test_case.name
        self.driver_signatures = driver_signatures
        self.run_dir = run_dir
        self.dependencies = test_dependencies
        self.ipi_command = shlex.split('i-pi' +
                                  ' ' + os.path.join(run_dir, 'input.xml'))

    def run(self):
        ipi_output_path = os.path.join(self.run_dir, Parameters.IPI_output_file)
        # Run the i-pi code
        os.chdir(self.run_dir)
        driver_prcs = []
        try:
            with open(ipi_output_path, 'w') as ipi_out:
                ipi_out.write('*REGTEST* IPI COMMAND: %s\n' % ' '.join(self.ipi_command))
                try:
                    ipi_proc = subprocess.Popen(self.ipi_command,
                                          bufsize=0,
                                          stdout=ipi_out,
                                          stderr=subprocess.PIPE)
                except OSError as _err:
                    if _err.errno == os.errno.ENOENT:
                        raise WrongIPICommand('i-pi command not found!', self.ipi_command)
                    else:
                        raise
            time.sleep(5.0)
            
            # Run the driver code
            for cmd, instance_folder in self.driver_signatures:
                os.chdir(instance_folder)
                driver_out_path = os.path.join(instance_folder, 'driver.out')
                with open(driver_out_path, 'w') as driver_out:
                    driver_out.write('*REGTEST* DRIVER COMMAND: %s\n' %
                                     ' '.join(cmd))
                    try:
                        #print cmd
                        driver_prcs.append(subprocess.Popen(cmd,
                                                      bufsize=0,
                                                      stdin=None,
                                                      stdout=driver_out,
                                                      stderr=subprocess.STDOUT))
                    except OSError as _err:
                        if _err.errno == os.errno.ENOENT:
                            raise WrongDriverCommand('driver command %s not found!\n' %
                                         ' ' .join(cmd))
                        else:
                            raise
                            
            driver_init_time = time.time()
            finished_drivers = 0

            while finished_drivers < len(driver_prcs):
                for prc in driver_prcs:
                    if prc.poll() is not None:
                        finished_drivers += 1
                time.sleep(Parameters.sleep_time)
                time_elapsed = time.time()-driver_init_time
                if time_elapsed > Parameters.driver_timeout:
                    for prc in driver_prcs:
                        if prc.poll() is None:
                            prc.terminate()
                    ipi_proc.terminate()
                    time.sleep(Parameters.ipi_shutdown_time)  # i-pi took approximately 10 sec to exit
                    stdout, stderr = ipi_proc.communicate()
                    if stderr:
                        with open(ipi_output_path, 'a') as ipi_out:
                            ipi_out.write(stderr)
                        raise IPIError("I-PI Error")
                    raise RegtestTimeout('Driver timeout')

            drivers_terminated_time = time.time()
            while ipi_proc.poll() is None:
                time.sleep(Parameters.sleep_time)
                time_elapsed_after_drivers = time.time() - drivers_terminated_time
                if time_elapsed_after_drivers > Parameters.ipi_shutdown_time:
                    ipi_proc.terminate()
                    time.sleep(Parameters.ipi_shutdown_time)  # i-pi took approximately 10 sec to exit
                    stdout, stderr = ipi_proc.communicate()
                    if stderr:
                        with open(ipi_output_path, 'a') as ipi_out:
                            ipi_out.write(stderr)
                        raise IPIError("I-PI Error")
                    raise RegtestTimeout('I-PI timeout after drivers finished')
                    
            stdout, stderr = ipi_proc.communicate()
            if stderr:
                with open(ipi_output_path, 'a') as ipi_out:
                    ipi_out.write(stderr)
                raise IPIError("I-PI Error")
        except KeyboardInterrupt:
            if ipi_proc.poll() is not None:
                ipi_proc.terminate()
                time.sleep(Parameters.ipi_shutdown_time)
                stdout, stderr = ipi_proc.communicate()
                if stderr:
                    with open(ipi_output_path, 'a') as ipi_out:
                        ipi_out.write(stderr)
            for prc in driver_prcs:
                    if prc.poll() is None:
                        prc.terminate()
            raise

    def get_output(self):
        return TestOutput(self.path_ipi_input, self.run_dir)
    def put_output_as_reference(self):
        output = self.get_output()
        output.put(self.test_case.reference_dir)
        return
    def compare_output_with_reference(self):
        output = self.get_output()
        reference = self.test_case.get_reference_output()
        return output.compare(reference)

class TestOutput:
    def __init__(self, xml, path):
        abs_xml = os.path.abspath(xml)
        if not os.path.isfile(abs_xml):
            raise ValueError("Input file does not exist %s" % abs_xml)
        #print "Before get_filesname"
        filelist = get_filesname(abs_xml)
        #print "After get_filesname"
        abs_files = [os.path.join(path, x) for x in filelist]
        filetuples = [(os.path.join(path, x),x) for x in filelist]
        for files in abs_files:
            if not (os.path.isfile(files)):
                raise OutputCorrupted("Expected output file %s not exist" % files)
        self.files = abs_files
        self.fileset = set(filelist)
        self.filetuples = set(filetuples)

    def compare(self, test_output):
        # Check if file lists are the same
        assert(self.fileset == test_output.fileset)
        # Check if they all exist
        for files in self.files:
            assert(os.path.isfile(files))
        for files in test_output.files:
            assert(os.path.isfile(files))
        #print zip(self.files,test_output.files)
        report = []
        for file1 in self.filetuples:
            file2 = [x for x in test_output.filetuples if x[1]==file1[1]]
            #print file1[0], file2[0][0]
        #for file1,file2 in zip(self.files,test_output.files):
            differences = compare_files(file1[0],file2[0][0])
            result = ComparisonResult(file1[0],file2[0][0], differences)
            report.append(result)
            #result.print_differences()
        return report
        
    def put(self, path):
        try:
            os.makedirs(path)
        except OSError:
            if os.path.exists(path):
                raise ValueError ("Directory %s exists. The dir parameter must be a path to nonexistent directory" % str(path))
            else:
                raise
        for file in self.files:
            shutil.copy(file, path)
        return

class ComparisonResult:
    def __init__(self, file1, file2, differences):
        self.file1 = file1
        self.file2 = file2
        self.list_of_differences = differences
    def append(self,place):
        self.list_of_differences.append(place)
    def print_differences(self):
        print "Files:", self.file1, self.file2
        for element in self.list_of_differences:
            print "Difference in line", element[0], "in word", element[1]
    def files_are_equal(self):
        return ( len(self.list_of_differences) == 0 )
        
class TestCandidate:
    def __init__(self, name, input_path):
        input_abspath = os.path.abspath(input_path)
        if file_is_test(os.path.join(input_abspath)):
            self.input_path = input_abspath
            self.name = name
        else:
            raise ValueError("Not a valid xml test")
            
class Parameters:
    regtest_string = r'<!--\s*REGTEST\s+([\s\w\W]*)\s+ENDREGTEST\s*-->'
    command_string = r'^COMMAND\(?(\d*)\)?\s*([ \w.\+\-\(\)\<\>\:]*)$'
    dependencies_string = r'^DEPENDENCIES\s*([ \w.\+\-\(\)]*)$'
    
    default_run_directory = 'regtest-run'
    default_test_cases_directory = '.'
    reference_directory = 'regtest-ref'
    ipi_timeout = 10
    default_driver_timeout = 600
    default_precision = 7
    precision = 7
    driver_timeout = 600
    IPI_output_file = 'ipi_output.out'
    ipi_shutdown_time = ipi.utils.softexit.SOFTEXITLATENCY * 1.5 # waits at least as much as the softexit latency before giving up hopes
    sleep_time = 0.5
    
    @staticmethod
    def compile_regtest_string():
        return re.compile(Parameters.regtest_string)
    @staticmethod
    def compile_command_string():
        return re.compile(Parameters.command_string, re.MULTILINE)
    @staticmethod
    def compile_dependencies_string():
        return re.compile(Parameters.dependencies_string, re.MULTILINE)
            
class RegtestTimeout(Exception):
    ''' Raise when subprocess of regttest timeout'''
    
class WrongIPICommand(Exception):
    ''' Raise when i-pi command is not found'''
    
class WrongDriverCommand(Exception):
    ''' Raise when driver command is not found'''
    
class OutputCorrupted(Exception):
    ''' Raise when test output (either reference or test) has missing files'''
    
class IPIError(Exception):
    ''' Raise when I-PI terminates with error'''

def inplace_change(filename, old_string, new_string):
    """ Replace a string in a file.

    Replace 'old_string' with 'new_string' into 'filename'.

    Args:
        filename: The filename where the string must be replaced.
        old_string: The string to replace.
        new_string: The string that will be replaced.

    Returns:
        Will be True if something has been replaced, False otherwise.
    """
    # Safely read the input filename using 'with'
    with open(filename) as _file:
        _content = _file.read()
        if old_string not in _content:
            return False
        if new_string in _content:
            return True

    # Safely write the changed content, if found in the file
    with open(filename, 'w') as _file:
        _content = _content.replace(old_string, new_string)
        _file.write(_content)
        return True

def get_filesname(xml_path):
        """
        Returns:
            lprop
            nprop
        """
        # Avoid to print i-pi output
        devnull = open('/dev/null', 'w')
        oldstdout_fno = os.dup(sys.stdout.fileno())
        os.dup2(devnull.fileno(), 1)
        
        xml_path = os.path.abspath(xml_path)
        os.chdir(os.path.dirname(xml_path))
        #I-pi xml file parser
        ifile = open(xml_path, "r")
        xmlrestart = io_xml.xml_parse_file(ifile)
        ifile.close()

        isimul = InputSimulation()
        isimul.parse(xmlrestart.fields[0][1])
        #print isimul.output
        simul = isimul.fetch()

        # reconstructs the list of the property and trajectory files
        lprop = []  # list of property files
        ltraj = []  # list of trajectory files
        for o in simul.outtemplate:
            # properties and trajectories are output per system
            if isinstance(o, CheckpointOutput):
                pass
            elif isinstance(o, PropertyOutput):
                for _ss in simul.syslist:   # create multiple copies
                    filename = _ss.prefix + o.filename
                    lprop.append(filename)

            # trajectories are more complex, as some have per-bead output
            elif isinstance(o, TrajectoryOutput):
                if getkey(o.what) in ["positions", "velocities",
                                      "forces", "forces_sc", "extras"]:   # multiple beads
                    nbeads = simul.syslist[0].beads.nbeads
                    for _bi in range(nbeads):
                        # zero-padded bead number
                        padb = (("%0" + str(int(1 +
                                                np.floor(np.log(nbeads) /
                                                         np.log(10)))) +
                                 "d") % (_bi))

                        for _ss in simul.syslist:
                            if o.ibead < 0 or o.ibead == _bi:
                                if getkey(o.what) == "extras":
                                    filename = _ss.prefix + o.filename + "_" + padb
                                else:
                                    filename = _ss.prefix + o.filename + "_" + padb + \
                                        "." + o.format
                                ltraj.append(filename)
                else:
                    for _ss in simul.syslist:   # create multiple copies
                        filename = _ss.prefix + o.filename
                        filename = filename + "." + o.format
                        ltraj.append(filename)
        os.dup2(oldstdout_fno, 1)
        return (ltraj + lprop)

def compare_files(file1,file2):
        """ Function to compare files.

        The strings are compared for equality and floats are compared
        to some precision using numpy.
        
        TODO: now integers are also converted to floats
        
        returns list of tuples of differences (line, word)
        """
        
        differences = []
        with open(file1, 'r') as content_file1:
            content1 = content_file1.readlines()
        with open(file2, 'r') as content_file2:
            content2 = content_file2.readlines()
        
        line_count = 1
        for line_in_file1, line_in_file2 in zip(content1, content2):
            word_count = 1
            for word_in_file1, word_in_file2 in zip(line_in_file1.split(),
                                    line_in_file2.split()):
                try:
                    float_in_file1 = float(word_in_file1)
                    float_in_file2 = float(word_in_file2)
                    if not np.isclose(float_in_file1, float_in_file2):
                        differences.append((line_count, word_count))
                        #print  float_in_file1, float_in_file2
                except ValueError:
                    if not word_in_file1 == word_in_file2:
                        #print word_in_file1, word_in_file2
                        differences.append((line_count, word_count)) 
                word_count += 1
            line_count += 1
        return differences

def answer_is_y():
    """ A simple function to interrogate the user on y/n questions.

    Only 'yes' and 'y' are counted as yes answer and only 'n' and 'no' are
    valied negative answer. All the answer are case insensitive. The function
    will ask for an answer until the user do not reply with a valid character.

    Args:
        msg: The question...

    Return:
        True if the user answer yes or False if the user answer no.
    """

    _yes = ['yes', 'y']
    _no = ['no', 'n']

    answer = raw_input()

    if answer.lower() in _yes:
        return True
    elif answer.lower() in _no:
        return False

if __name__ == '__main__':
    main()
