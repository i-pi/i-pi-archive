.PHONY: all clean high_density mid_density low_density
all: high_density mid_density low_density

driver:=../../../drivers/driver.x
IPI:=../../../i-pi

define run_driver
  for i in `seq 1 $1`; do \
    $(driver) -u -h $2 -m lj -o 5.270446,1.1663e-4,13.176115 & \
  done;
endef
-include make.in

high_density:
	cd high_density; $(IPI) input.xml & sleep 5; \
	$(call run_driver,1,high_density) \
	wait
	
mid_density:
	cd mid_density; $(IPI) input.xml & sleep 5; \
	$(call run_driver,1,mid_density) \
	wait

low_density:
	cd low_density; $(IPI) input.xml & sleep 5; \
	$(call run_driver,1,low_density) \
	wait

clean:
	rm -f */lj-* */RESTART */EXIT
