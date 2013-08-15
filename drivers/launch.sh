#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe  orte 4
#$ -N lj-ipi

# Script to run an individual RPMD run
# 
# Copyright (C) 2013, Joshua More and Michele Ceriotti
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http.//www.gnu.org/licenses/>.


port=$1
nprocs=4

for i in `seq 1 $nprocs`; do
   ./driver.x -m sg -o 15.0 -h localhost -p $port &
done
wait
