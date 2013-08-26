# Script to run one of the directories of RPMD tests
# 
# Copyright (C) 2013, Joshua More and Michele Ceriotti
# 
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
# 
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

nruns=99
port=$1

python ../../../../i-pi input.xml &> log1
wait
for i in `seq 1 $nruns`; do
   newport=$((5+$port))
   j=$(($i+1))
   sed "s/<simulation>/<simulation> <initialize nbeads='24'> <velocities mode='thermal' units='kelvin'> 25 <%velocities> <%initialize>/" RESTART | tr "%" "/" > dummy
   sed "s/<step>8000</<step>0</; s/<port>$port/<port>$newport/; s/filename='test$i/filename='test$j/" dummy > input$i

   port=$(($newport)) 

   rm dummy
   python ../../../../i-pi input$i &> log$j
   wait
done
