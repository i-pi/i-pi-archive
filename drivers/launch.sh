#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe  orte 4
#$ -N lj-wrappi

port=$1
nprocs=4

for i in `seq 1 $nprocs`; do
   ./driver.x -m sg -o 15.0 -h localhost -p $port &
done
wait
