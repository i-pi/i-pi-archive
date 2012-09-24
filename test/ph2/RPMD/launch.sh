nruns=99
port=$1

python ../../../../src/main.py input.xml &> log1
wait
for i in `seq 1 $nruns`; do
   newport=$((5+$port))
   j=$(($i+1))
   sed "s/<simulation>/<simulation> <initialize> <rescale_v> 25 <%rescale_v> <%initialize>/" test$i.restart1 | tr "%" "/" > dummy
   sed "s/<step>8000</<step>0</; s/<port>$port/<port>$newport/; s/prefix='test$i/prefix='test$j/" dummy > input$i

   port=$(($newport)) 

   rm dummy
   python ../../../../../src/main.py input$i &> log$j
   wait
done
