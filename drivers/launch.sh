#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe  orte 7
#$ -N lj-wrappi

port=$1
nprocs=6

for i in `seq 1 $nprocs`; do
   ./driver.x -m sg -o 15.0 -h neptune.chem.ox.ac.uk -p $port &
done
wait
