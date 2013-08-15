.PHONY: all clean npt-classical
all: npt-classical

CP2K:=cp2k.x
IPI:=../../../i-pi

define run_cp2k
  for i in `seq 1 $1`; do \
    mkdir -p run_$$i; cd run_$$i; \
    $(CP2K) -i ../$2 -o cp2k.out & cd ..; \
  done;
endef
-include make.in

npt-classical:
	cd npt-classical; $(IPI) nptcls_ipi.xml & sleep 3; \
	$(call run_cp2k,1,nptcls_cp2k.in) \
	wait
	
clean:
	rm -rf */*_output.* */*.log */run_* */RESTART */EXIT
