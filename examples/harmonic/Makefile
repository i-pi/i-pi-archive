.PHONY: all clean harmonic
all: harmonic

driver:=../../drivers/driver.x
IPI:=../../i-pi

define run_driver
  for i in `seq 1 $1`; do \
    $(driver) -u -h harmonic -m harm -o 1 & \
  done;
endef
-include make.in

harmonic:
	$(IPI) input.xml & sleep 5; \
	$(call run_driver,4) \
	wait

clean:
	rm -f harm.* RESTART EXIT
