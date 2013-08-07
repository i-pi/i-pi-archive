# Script to run one of the directories of RPMD tests
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
port=$1

python ../../../../src/main.py input.xml &> log1
wait
for i in `seq 1 $nruns`; do
   newport=$((5+$port))
   j=$(($i+1))
   sed "s/<simulation>/<simulation> <initialize> <resample_v> 25 <%resample_v> <%initialize>/" RESTART | tr "%" "/" > dummy
   sed "s/<step>8000</<step>0</; s/<port>$port/<port>$newport/; s/filename='test$i/filename='test$j/" dummy > input$i

   port=$(($newport)) 

   rm dummy
   python ../../../../src/main.py input$i &> log$j
   wait
done
