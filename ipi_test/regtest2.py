#!/usr/bin/env python2

import argparse
import xml.etree.ElementTree as etree
import os
import Queue
import re
import shutil
import shlex
import subprocess as sbps
import sys
import time
import threading

import numpy as np
import numpy.testing as npt

# pylint: disable=import-error
from ipi.engine.outputs import CheckpointOutput, PropertyOutput, \
    TrajectoryOutput
from ipi.engine.properties import getkey
from ipi.inputs.simulation import InputSimulation
from ipi.utils.io.inputs import io_xml
# pylint: enable=import-error


# Compile them only once! pylint: disable=anomalous-backslash-in-string
DRIVER_COMMAND_RGX = re.compile(r'\<\!\-\-\s*[Rr]un\s*the\s*driver\s*code'
                                '\s*with\:\s*"(.*)"\s*\-\-\>', re.MULTILINE)
NEEDED_FILES_RGX = re.compile(r'\<\!\-\-\s*[Tt]he\s*driver\s*will\s*need\s*the'
                              '\s*following\s*files\:\s*(.*)\s*\-\-\>',
                              re.MULTILINE)
# pylint: enable=anomalous-backslash-in-string

QUEUE_ALL = Queue.Queue() # Queue to access the "run"
QUEUE_COM = Queue.Queue() # Queue to access the "compare"



def main():

    root_test_folder = _parser()['root_test_folder']
    tests_list = _parser()['tests']

    # Create the root folder of the run
    root_run = _parser()['root_run_folder']
    create_dir(root_run)


    # If no --add-test has been specified, search for tests in all directories
    #+within the root_test_folder.
    if len(tests_list) == 0:
        tests_list = ['']

    # Retrieve the test list and build the QUEUE_ALL
    test_list = _build_test_index(root_test_folder, tests_list)
    index = _parser()['index']
    for test in test_list:
        test_obj = Test(index = index, name=test[0], path=test[1], root_run=root_run)
        test_obj.run()
        # test_obj.create_reference()
        test_obj.run()
        test_obj.print_report()
    #     QUEUE_ALL.put(test_obj)
    #     index += 10
    # running_test = []
    # while True:
    #     if len(running_test) < _parser()['nproc']:
    #         try:
    #             running_test.append(QUEUE_ALL.get_nowait())
    #         except Queue.Empty:
    #             pass
    #         else:
    #             running_test[-1].start()

    #     for _ii, thr in enumerate(running_test):
    #         thr.join(0.5)
    #         if thr.is_alive():
    #             QUEUE_COM.put(thr)
    #             QUEUE_ALL.task_done()
    #             del running_test[_ii]

    #     if not _parser()['create_reference']:
    #         try:
    #             thr_comparing = QUEUE_COM.get_nowait()
    #         except Queue.Empty:
    #             pass
    #         else:
    #             thr_comparing.start()
    #             thr_comparing.join(0.5)
    #             if thr_comparing.is_alive():
    #                 thr_comparing.print_report()
    #                 QUEUE_COM.task_done()

    #     if len(running_test) == 0 and QUEUE_ALL.empty() and QUEUE_COM.empty():
    #         break



def _parser():
    """ Parse the argument lists given as input.

    Return:
        A dictionary containing all the input options and argument.
    """
    parser = argparse.ArgumentParser(description='Regtests for i-PI.')
    parser.add_argument('--add-test',
                        action='append',
                        type=str,
                        help=('Mark a test to be ran. A single test can be '
                              'specified with "group:test". This option can be '
                              'used as many times as you want.'),
                        default=[],
                        dest='tests')
    parser.add_argument('--tests-folder',
                        action='store',
                        type=str,
                        default='example',
                        help=('The folder where to search for tests.'),
                        dest='root_test_folder')
    parser.add_argument('--folder-run',
                        action='store',
                        type=str,
                        default='regtest-run',
                        help=('Parent folder where the test will be run. '
                              'If already existing all the content will '
                              'be lost!'),
                        dest='root_run_folder')
    parser.add_argument('--starting-address',
                        action='store',
                        type=int,
                        default=10000,
                        help=('Number of the first port in case of inet '
                              'socket.'),
                        dest='index')
    parser.add_argument('-np', '--nproc',
                        action='store',
                        type=int,
                        default=4,
                        help=('Number of concurrent test run at once.'),
                        dest='nproc')
    parser.add_argument('--create-reference',
                        action='store_true',
                        default=False,
                        help=('Only run the test, do not compare and at the '
                              'end ask if you want to copy the output files '
                              'in the reference test folder.'),
                        dest='create_reference')
    parser.add_argument('--precision',
                        action='store',
                        default=7,
                        type=int,
                        help=('Define the number of decimal used in comparing'
                              'float.'),
                        dest='precision')


    return vars(parser.parse_args())



def _build_test_index(root_path, tests):
    """ Look for a valid xml in all the specified directories.

    Check if there are valid xml in all the specified folders. In case append
    the test to a list. If a folder contains more than a single xml, that folder
    will be skipped and a warning is printed. If a folder do not contain any xml
    that folder will be skipped without any message.
    Once all the folder are examinated, a list of all the found tests is printed
    at stdout.

    Print to stdout a list of all the tests found.

    Args:
        root_path: Parent folder of all the tests.
        tests: A sequence of subfolder of root_path to filter the number of
            tests found.

    Returns:
        test_list: A sequence containing all the valid tests found. Each
            element is a tuple of the form (<test_name>, <abs_path_to_xml>)

    """

    root_path = os.path.abspath(root_path)
    tests_folder = [os.path.join(root_path, test) for test in tests]
    test_list = []
    msg = '**Tests that will be executed:**\n'

    for test in tests_folder:
        if not os.path.exists(test) or not os.path.isdir(test):
            raise RuntimeError('Folder %s does not exist!' % test)
        for root, junk, files in os.walk(test): # pylint: disable=unused-variable
            test_name = os.path.relpath(root, root_path)
            xml_files = [x for x in files if x.endswith('.xml')]
            if len(xml_files) > 1:
                sys.stderr.write('!W! Skipping test %s: too many xml files!\n'\
                                 % test_name)
                continue
            elif len(xml_files) < 1:
                continue

            if _file_is_test(os.path.join(root, xml_files[0])):
                test_list.append((test_name,
                                  os.path.join(root, xml_files[0])))
                msg += ' > ' + str(os.path.split(root)[1]) + '\n'

    if len(test_list) < 1:
        print "**No test found!**"
    else:
        print msg

    return test_list




def _file_is_test(path_to_test):
    """ Check if an xml file can be used as regtest.

    To be a valid regtest the xml file must contain the line to specify which
    is the command to be executed as a driver.
    """

    with open(path_to_test) as _file:
        _text = _file.read()
    return len([x.group(1) for x in DRIVER_COMMAND_RGX.finditer(_text)]) > 0


#class Test(threading.Thread):
class Test(object):

    def __init__(self, *args, **kwds):
#        threading.Thread.__init__(self)

        self.index = kwds['index']
        self.test_path = kwds['path']
        self.name = kwds['name']
        self.run_path = os.path.abspath(kwds['root_run'])

        del kwds['index']
        del kwds['path']
        del kwds['name']
        del kwds['root_run']

        self.ref_path = os.path.dirname(self.test_path)
        self.filename = os.path.basename(self.test_path)
        self.msg = ''
        self.driver_command = None
        self.input_dir = None
        self.io_dir = None
        self._test_status = 'PASSED'

        self.generate_output = True
        self.compare_output = False

    @property
    def test_status(self):
        """ Status of the test: always print the most important status!

        An error has always the precedence, then a failure and, only if nothing
        bad happen, print PASSED.
        """
        return self._test_status

    @test_status.setter
    def test_status(self, status):
        status_priority = {'ERROR': 10,
                           'FAILED': 5,
                           'PASSED': 1}

        if status_priority[status] > status_priority[self._test_status]:
            self._test_status = status


    def init_env(self):
        """ Prepare the test run.

        Do the following thing:
        * Create the necessary folders.
        * Retrieve the command to run the driver from the xml.
        * Copy all the reference files into the 'input' folder.
        * Copy all the necessary files to the 'io' folder.
        * Change the address of the socket to be unique

        There is some rule on the way the socket address is changed and it is
        important to know this rules when preaparing the test:
        * If unix socket then the <address> tag is changed appending a number
          defined using the self.index.
        * If inet socket then the <port> tag to a number defined using the
          self.index.
        * To ensure that the driver is using the right address the first word
          of the driver command that match the address of the xml file is
          replaced with the new address.
        """

        # Create the necessary folders
        path = self.run_path
        for folder in self.name.split('/'):
            path = os.path.join(path, folder)
            create_dir(path, ignore=True)
        self.run_path = path

        self.io_dir = os.path.join(self.run_path, 'io')
        self.input_dir = os.path.join(self.run_path, 'input')
        create_dir(self.io_dir, ignore=True)

        # Retrieve the command to run the driver and the needed files
        with open(self.test_path) as _buffer:
            _text = _buffer.read()

        self.driver_command = [x.group(1) for x in
                               DRIVER_COMMAND_RGX.finditer(_text)]
        _buffer = [x.group(1) for x in NEEDED_FILES_RGX.finditer(_text)]
        needed_files = [os.path.basename(self.test_path)]
        for _xx in _buffer:
            needed_files += [x.strip() for x in _xx.split(',')]

        # Copy all the reference files to the 'input' folder
        shutil.copytree(os.path.dirname(self.test_path), self.input_dir)

        # Copy all only the needed files to the 'io' folder
        for _buffer in needed_files:
            shutil.copy2(os.path.join(self.input_dir, _buffer), self.io_dir)

        # Make sure to do not change the reference
        self.test_path = os.path.join(self.io_dir, self.filename)

        # Change the address of the socket to be unique
        xml = etree.parse(self.test_path)
        address_index = self.index
        for ffsocket in xml.findall('./ffsocket'):
            socket_mode = ffsocket.attrib['mode']
            if socket_mode.lower() == 'unix':
                address = ffsocket.find('./address').text.strip()
                ffsocket.find('./address').text = (' %s%i ' % (address,
                                                               address_index))

                # Change the address used by the driver too!
                for _ii, _buffer in enumerate(self.driver_command):
                    # Replace only the first occurrence of 'address' in the
                    #+driver command!
                    self.driver_command[_ii] = _buffer.replace(address,
                                                               ' %s%i ' % \
                                                               (address,
                                                                address_index),
                                                               1)
            # Do the same when using inet socket
            else:
                address = ffsocket.find('./port').text.strip()
                ffsocket.find('./port').text = str(address_index)
                for _ii, _buffer in enumerate(self.driver_command):
                    self.driver_command[_ii] = _buffer.replace(address,
                                                               ' %i ' % \
                                                               (address_index),
                                                               1)

            # Since there could be more than a single socket used within a
            #+single simulation, it is important to have a different address for
            #+each socket.
            address_index += 1

        xml.write(self.test_path)

        # Change to the io directory...
        os.chdir(self.io_dir)


    def run(self):
        if self.generate_output:
            self.init_env()
            self._run_ipi()
            self.generate_output = False
            self.compare_output = True
        elif self.compare_output:
            self._compare()
            self.compare_output = False


    def _run_ipi(self):
        # Move to the right directory

        driver_out_path = os.path.join(self.io_dir, 'driver_output.out')

        ipi_command = 'i-pi'
        timeout_driver = 5
        timeout_ipi = 5

        # Run the i-pi code
        ipi_command = shlex.split(ipi_command +\
                                  ' ' + self.test_path)
        with open(os.path.join(self.io_dir, 'ipi_output.out'), 'w') as ipi_out:
            try:
                ipi_proc = sbps.Popen(ipi_command,
                                      bufsize=0,
                                      stdout=ipi_out,
                                      stderr=sbps.PIPE)
            except OSError as _err:
                if _err.errno == os.errno.ENOENT:
                    raise OSError('i-pi command not found!')

        # Sleep few seconds waiting for the ipi server start
        time.sleep(2)

        # Run the driver code
        driver_prcs = []
        for cmd in self.driver_command:
            cmd = shlex.split(cmd)
            with open(driver_out_path, 'w') as driver_out:
                try:
                    driver_prcs.append(sbps.Popen(cmd,
                                                  bufsize=0,
                                                  stdin=None,
                                                  stdout=driver_out,
                                                  stderr=sbps.STDOUT))
                except OSError as _err:
                    if _err.errno == os.errno.ENOENT:
                        raise OSError('driver command %s not found!' %
                                      ' ' .join(cmd))

        init_time = -time.time()
        finished = 0
        while finished != len(driver_prcs):
            for prc in driver_prcs:
                if prc.poll() is not None:
                    finished += 1
            time.sleep(.5)
            timeout_driver = timeout_driver - init_time - time.time()
            if timeout_driver < -5:
                for prc in driver_prcs:
                    prc.kill()
                    finished += 1
                self.test_status = 'ERROR'
                self.msg += '\nThe drivers took too long!'

        while ipi_proc.poll() is None:
            timeout_ipi -= 2
            time.sleep(2)
            if timeout_ipi < -2:
                ipi_proc.kill()
                self.test_status = 'ERROR'
                self.msg += '\ni-PI took too long after the driver finished!'

        stdout, stderr = ipi_proc.communicate() # pylint: disable=unused-variable
        if len(stderr) != 0:
            self.test_status = 'ERROR'
            self.msg += '\ni-PI produced the following error:\n'
            self.msg += stderr
            self.msg += '\n'

        return

    def create_reference(self):

        if self.test_status == 'PASSED':
            if answer_is_y('Do you want to use the produced files as reference and copy them on the original folder?\nThis will overwrite any content of the original folder with the test/io folder!'):

                self.input_dir = self.ref_path
                ltraj, lprop =  self._get_filesname()

                for prop in lprop:
                    for sprop in prop:
                        remove_file(sprop['old_filename'])
                        shutil.copy2(sprop['new_filename'], sprop['old_filename'])
                for traj in ltraj:
                    for straj in traj:
                        remove_file(straj['old_filename'])
                        shutil.copy2(straj['new_filename'], straj['old_filename'])



    def print_report(self):
        _format = '%30s -->  %15s Info: %s\n'
        msg = _format % (self.name, self.test_status, self.msg)
        print msg


    def _get_filesname(self):
        """ The test results shold be analyzed number by numbers.

        The idea is that the testing input should never change, then the files
        should be always the same. It would probably be better, anyway, to use the
        i-pi infrastructure to retrieve the right position of the data. In fact,
        this would work as a further testing.

        Returns:
            lprop
            nprop
        """

        # opens & parses the input file
        ifile = open(self.test_path, "r")
        xmlrestart = io_xml.xml_parse_file(ifile) # Parses the file.
        ifile.close()

        isimul = InputSimulation()
        isimul.parse(xmlrestart.fields[0][1])

        simul = isimul.fetch()

        # reconstructs the list of the property and trajectory files
        lprop = [] # list of property files
        ltraj = [] # list of trajectory files
        for o in simul.outtemplate:
            # properties and trajectories are output per system
            if isinstance(o, CheckpointOutput):
                pass
            elif isinstance(o, PropertyOutput):
                nprop = []
                isys = 0
                for _ in simul.syslist:   # create multiple copies
                    filename = o.filename
                    nprop.append({"old_filename" : os.path.join(self.input_dir, filename),
                                  "new_filename" : os.path.join(self.io_dir, filename),
                                  "stride": o.stride,
                                  "properties": o.outlist,})
                    isys += 1
                lprop.append(nprop)

            # trajectories are more complex, as some have per-bead output
            elif isinstance(o, TrajectoryOutput):
                if getkey(o.what) in ["positions", "velocities",
                                      "forces", "extras"]:   # multiple beads
                    nbeads = simul.syslist[0].beads.nbeads
                    for _bi in range(nbeads):
                        ntraj = []
                        isys = 0
                        # zero-padded bead number
                        padb = (("%0" + str(int(1 +
                                                np.floor(np.log(nbeads) /
                                                         np.log(10)))) +
                                 "d") % (_bi))

                        for _ in simul.syslist:
                            if o.ibead < 0 or o.ibead == _bi:
                                if getkey(o.what) == "extras":
                                    filename = o.filename+"_" + padb
                                else:
                                    filename = o.filename+"_" + padb + \
                                               "." + o.format
                                ntraj.append({"old_filename" : os.path.join(self.input_dir, filename),
                                              "format" : o.format,
                                              "new_filename" : os.path.join(self.io_dir, filename),
                                              "stride": o.stride,
                                              "what": o.what,})
                            isys += 1
                        if ntraj != []:
                            ltraj.append(ntraj)

                else:
                    ntraj = []
                    isys = 0
                    for _ in simul.syslist:   # create multiple copies
                        filename = o.filename
                        filename = filename+"."+o.format
                        ntraj.append({"old_filename" : os.path.join(self.input_dir,
                                                                    filename),
                                      "new_filename" : os.path.join(self.io_dir,
                                                                    filename),
                                      "format" : o.format,
                                      "stride": o.stride,})

                        isys += 1
                    ltraj.append(ntraj)

        return ltraj, lprop


    def compare_property_files(self, lprop):
        """ Comparing property files.

        The files are loaded with np.loadtxt and all the intestation are ignored:
        only the actual values are compared. Thus, the ipi input must define the
        same columns in the same order.
        """

        for prop in lprop:
            for sprop in prop:
                try:
                    old_content = np.loadtxt(sprop['old_filename'])
                except IOError as _err:
                    if _err.errno == os.errno.ENOENT:
                        self.msg += 'File %s not found!\n' % sprop['old_filename']
                        self.test_status = 'ERROR'
                        continue

                try:
                    new_content = np.loadtxt(sprop['new_filename'])
                except IOError as _err:
                    if _err.errno == os.errno.ENOENT:
                        self.msg += 'File %s not found!\n' % sprop['old_filename']
                        self.test_status = 'ERROR'
                        continue

                try:
                    npt.assert_array_almost_equal(old_content, new_content,
                                                  _parser()['precision'])
                except AssertionError:
                    name = os.path.basename(sprop['old_filename'])
                    self.msg +=  'Differences in the %s file\n' % name
                    self.test_status = 'FAILED'
                    continue



    def compare_trajectory_files(self, ltraj):
        """ Function to compare trajectory files.

        The idea is to store all the numbers in the file in a list and then using
        numpy to compare the two lists. The numbers are recognized exploiting the
        float function error when applied on strings.

        The strings are compared directly.
        """

        err = False
        for traj in ltraj:
            for straj in traj:
                new_w_list = []
                old_w_list = []
                name = os.path.basename(straj['old_filename'])
                try:
                    old_content = open(straj['old_filename'])
                except IOError as _err:
                    if _err.errno == os.errno.ENOENT:
                        self.msg += 'File %s not found!\n' % straj['old_filename']
                        self.test_status = 'ERROR'
                        continue

                try:
                    new_content = open(straj['new_filename'])
                except IOError as _err:
                    if _err.errno == os.errno.ENOENT:
                        self.msg += 'File %s not found!\n' % straj['new_filename']
                        self.test_status = 'ERROR'
                        continue


                line_c = 1
                for old_line, new_line in zip(old_content, new_content):
                    word_c = 1
                    for old_w, new_w in zip(old_line.split(),
                                            new_line.split()):
                        try:
                            old_w_list.append(float(old_w))
                            new_w_list.append(float(new_w))
                        except ValueError:
                            try:
                                assert old_w == new_w
                            except AssertionError:
                                self.msg += 'Differences at line %d word %d of file  %s' % (line_c, word_c, name)
                                self.test_status = 'FAILED'

                        word_c += 1
                    line_c += 1

                try:
                    npt.assert_array_almost_equal(np.array(new_w_list),
                                                  np.array(old_w_list),
                                                  _parser()['precision'])
                except AssertionError:
                    self.msg += 'Differences in the %s file\n' % name
                    self.test_status = 'FAILED'
                    continue

        return err



    def _compare(self):
        """ This is the function that compares all the ipi output.

        The name of the files to compare come from ltraj and lprop.
        """

        ltraj, lprop =  self._get_filesname()

        self.compare_property_files(lprop)

        self.compare_trajectory_files(ltraj)




#######################
### Tools Functions ###
#######################

def remove_file(path):
    """ Remove a path only if it exists and is a file!

    """
    if os.path.exists(path) and os.path.isfile(path):
        os.remove(path)

def create_dir(folder_path, ignore=False):
    """ Create a folder after user consense.

    If the path asked by the user is not existing, just create that folder
    otherwise asks the user his opinion on deleting the existing folder/file.

    Args:
        folder_path: The path of the folder that should be created.
        ignore: If True and the path already exists, go on...

    """
    if os.path.exists(folder_path):
        if ignore:
            return True
        if answer_is_y('!W! %s already exists!\nDo you want to delete it and '
                       'proceed?\n' % folder_path):
            try:
                if os.path.isdir(folder_path):
                    shutil.rmtree(folder_path)
                elif os.path.isfile(folder_path):
                    os.remove(folder_path)
                else:
                    raise RuntimeError
            except:
                raise RuntimeError('I cannot remove the file.'
                                   'Try manually and restart this script!')
        else:
            raise SystemExit('User rules!')

    os.mkdir(folder_path)

    return True


def answer_is_y(msg):
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
    if not msg.endswith('\n'):
        msg += '\n'

    answer = ''
    while answer.lower() not in _yes and answer not in _no:
        sys.stderr.write(msg)

        answer = raw_input()

        if answer.lower() in _yes:
            return True
        elif answer.lower() in _no:
            return False



if __name__ == '__main__':
    main()
