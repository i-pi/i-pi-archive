Using driver.x with (high quality) potential energy surfaces
============================================================

Zundel cation
-------------

The `zundel/` folder contains a simple example of how to run driver.x
together with a potential fit for the Zundel cation in the gas
phase based on high-end quantum chemistry methods. 

One simply needs to run i-pi with

```bash
../../../i-pi input.xml
```

and then run one or more instances of the driver code, that 
contains a routine by Xinchuan Huang, Bastiaan J. Braams, and 
Joel M. Bowman, [J. Chem. Phys. 122, 044308 (2005)] to compute 
the energy for H5O2+. Note that the data files `h5o2.pes4B.coeff.dat`
and `h5o2.dms4B.coeff.com.dat` should be present in the folder
one wants to run the driver code. 

```bash
../../../drivers/driver.x -u -h zundel -m zundel
```


qTIP4P/f water model
--------------------

The `qtip4pf/` folder contains an example of how to run driver.x
with a simplistic implementation of the qtip4pf/ water model. 

One as usual runs i-PI, followed by one or more instances of the driver:

```bash
../../../i-pi input.xml &
../../../drivers/driver.x -u -h driver -m qtip4pf
```
