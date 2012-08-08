nruns=99
for i in `seq 1 5`; do 
   cd run_$i
   rm log* test* EXIT
   for j in `seq 1 $nruns`; do 
      rm input$j
   done
   cd .. 
done
