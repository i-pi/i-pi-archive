.PHONY: all clean nvt
all: nvt

driver:=../../../drivers/driver.x
IPI:=../../../i-pi

nvt:
	$(IPI) input.xml & sleep 5; \
	$(driver) -u -h ph2_nvt -m sg -o 15.0; \
	wait

clean:
	rm ph2.* RESTART EXIT
