i-PI
--------

A Python interface for ab initio path integral molecular dynamics simulations. 



Quick Installation and Test 
---------------------------

Follow these instruction to test i-PI. These assume to be run from a Linux 
environment, with a recent version of Python, Numpy and gfortran, and that 
the terminal is initially in the i-pi package directory (the directory containing
this file).

* Generate the driver code

$ cd driver
$ make
$ cd ..

* Run one of the examples

This will first start the wrapper in the background, redirecting the output on 
a log file, then run a couple of instances of the driver code and then follow
the progress of the wrapper by monitoring the log file

$ cd examples/tutorial/tutorial-1/
$ ../../../i-pi tutorial-1.xml > log &
$ ../../../drivers/driver.x -h localhost -p 31415 -m sg -o 15 &
$ ../../../drivers/driver.x -h localhost -p 31415 -m sg -o 15 &
$ tail -f log

The monitoring can be interrupted with CTRL+C when the run has finished (5000 steps)

