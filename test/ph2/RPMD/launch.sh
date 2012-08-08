nruns=99
port=$1

python ../../../../src/main.py input.xml &> log1
wait
for i in `seq 1 $nruns`; do
   newport=$((5+$port))
   j=$(($i+1))
   sed "s/<simulation>/<simulation> <initialize> {velocities : 25} <%initialize>/g" test$i.restart1 | tr "%" "/" > dummy
   sed "s/<step>8000</<step>0</g; s/<port>$port/<port>$newport/g; s/<prefix>test$i/<prefix>test$j/g" dummy > input$i

   port=$(($newport)) 

   rm dummy
   python ../../../../../src/main.py input$i &> log$j
   wait
done
