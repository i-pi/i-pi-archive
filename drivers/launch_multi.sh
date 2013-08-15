#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe  orte 4
#$ -N ph2-ipi

# Script to run the RPMD test cases

nruns=100
port=$1
nprocs=4

for i in `seq 1 $nruns`; do
   bash ./launch.sh $port
   wait
   port=$((5+$port))
   sleep 10
done
