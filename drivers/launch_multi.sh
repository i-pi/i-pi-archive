#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe  orte 7
#$ -N ph2-wrappi

nruns=100
port=$1
nprocs=6

for i in `seq 1 $nruns`; do
   bash ./launch.sh $port
   wait
   port=$((5+$port))
   sleep 10
done
