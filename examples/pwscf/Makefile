.PHONY: all water lithium clean
all: water lithium

PW:=pw-driver.x
IPI:=../../i-pi

define run_qe
  for i in `seq 1 $1`; do \
    mkdir -p run_$$i; cd run_$$i; \
    $(PW) < ../$2 > qe.out & cd ..; \
  done;
endef
-include make.in

water:
	$(IPI) h2o-ipi.xml & sleep 3; \
	$(call run_qe,4,h2o-pw.in) \
	wait

lithium:
	$(IPI) li4-ipi.xml & sleep 3; \
	$(call run_qe,4,li4-pw.in) \
	wait
	
clean:
	rm -rf h2o_output.* li4_output.* *.log RESTART EXIT run_*
