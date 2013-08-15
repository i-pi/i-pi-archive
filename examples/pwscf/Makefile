.PHONY: all water lithium clean
all: water lithium

PW:=pw.x
IPI:=../../../i-pi

define run_qe
  for i in `seq 1 $1`; do \
    mkdir -p run_$$i; cd run_$$i; \
    $(PW) < ../$2 > qe.out & cd ..; \
  done;
endef
-include make.in

water:
	cd h2o; $(IPI) h2o-ipi.xml & sleep 3; \
	$(call run_qe,4,h2o-pw.in) \
	wait

lithium:
	cd li4; $(IPI) li4-ipi.xml & sleep 3; \
	$(call run_qe,4,li4-pw.in) \
	wait
	
clean:
	rm -rf */*_output.* */*.log */run_* */RESTART */EXIT
