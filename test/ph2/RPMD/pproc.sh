#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe  orte 7
#$ -N ph2-pproc

nruns=100

for j in `seq 1 5`; do
   cd run_$j
   for i in `seq 1 $nruns`; do
      ../vel_est $i &
   done
   wait
   cd ..
done
wait
./consolidate &
