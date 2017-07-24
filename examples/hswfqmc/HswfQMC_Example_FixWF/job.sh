#!/bin/bash
 
#SBATCH -J HswfQMC-PIMD_Example  # Job name
#SBATCH -o %j.out                # Specify stdout output file
#SBATCH -p nodeshort             # Queue name
#SBATCH -N 1                     # Total number of nodes requested (64 cores/node)
#SBATCH -n 64                    # Total number of tasks
#SBATCH -t 5:00:00               # Run time (hh:mm:ss)
 
#SBATCH -A hydroqmd              # Specify allocation to charge against
 
# Launch the executable
HswfQMC_ipic -h 89.163.249.253 -p 54321 -m fix -o 0 -c 64,1 -l 10 -v
