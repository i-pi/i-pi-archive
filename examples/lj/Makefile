.PHONY: all clean high_density mid_density low_density
all: high_density mid_density low_density

driver:=../../../drivers/driver.x
IPI:=../../../i-pi

high_density:
	cd high_density; $(IPI) input.xml & sleep 5; \
	$(driver) -u -h high_density -m lj -o 5.270446,1.1663e-4,13.176115; \
	wait
	
mid_density:
	cd mid_density; $(IPI) input.xml & sleep 5; \
	$(driver) -u -h mid_density -m lj -o 5.270446,1.1663e-4,13.176115; \
	wait

low_density:
	cd low_density; $(IPI) input.xml & sleep 5; \
	$(driver) -u -h low_density -m lj -o 5.270446,1.1663e-4,13.176115; \
	wait

clean:
	rm */lj-* */RESTART */EXIT
