# Script to remove all the files created in the course of the 
# RPMD simulations
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

nruns=99
for i in `seq 1 5`; do 
   cd run_$i
   rm -f log* test* EXIT RESTART
   for j in `seq 1 $nruns`; do 
      rm -f input$j
   done
   cd .. 
done
