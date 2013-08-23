.PHONY: all clean harmonic
all: harmonic

driver:=../../drivers/driver.x
IPI:=../../i-pi

harmonic:
	$(IPI) input.xml & sleep 5; \
	$(driver) -u -h harmonic -m harm -o 1; \
	wait

clean:
	rm -f harm.* RESTART EXIT
