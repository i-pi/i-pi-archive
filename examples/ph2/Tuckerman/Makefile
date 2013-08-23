.PHONY: all clean 25K 20K 17K 14K
all: 25K 20K 17K 14K

driver:=../../../../drivers/driver.x
IPI:=../../../../i-pi

25K:
	cd 25K; $(IPI) input_equilibrated.xml & sleep 5; \
	$(driver) -u -h T_25K -m sg -o 15.0; \
	wait
	
20K:
	cd 20K; $(IPI) input_equilibrated.xml & sleep 5; \
	$(driver) -u -h T_20K -m sg -o 15.0; \
	wait

17K:
	cd 17K; $(IPI) input_equilibrated.xml & sleep 5; \
	$(driver) -u -h T_17K -m sg -o 15.0; \
	wait

14K:
	cd 14K; $(IPI) input_equilibrated.xml & sleep 5; \
	$(driver) -u -h T_14K -m sg -o 15.0; \
	wait

clean:
	rm -f */ph2-* */RESTART */EXIT
