FC=pgfortran
FOPT= -O3
NUMDIR = ./numerical
QAGDIR = ./numerical/dqag
WDIR = ./Wfunctions
RDIR = ./Rfunctions
RDIRDEV = ./Rfunctions_dev

MAIN = mainOper.o
MFSHR = sharedcap.o
MFOBJ = gencap.o
MFGPU = gpucap.o
MFCAP = opercap.o
NUMOBJ =  dgamic.o d1mach.o
QAG=  dsntdqagse.o dqelg.o dqk21.o dqpsrt.o dsntdqk21.o
WFUNC = WM.o WS2.o WS1.o WP2.o WMP2.o WP1.o WD.o WS1D.o
RFUNC = RM.o RS2.o RS1.o RP2.o RMP2.o RP1.o RD.o RS1D.o
RFUNCDEV = RM_dev.o RS2_dev.o RS1_dev.o RP2_dev.o RMP2_dev.o RP1_dev.o RD_dev.o RS1D_dev.o


gentest.x: $(MAIN) $(MFSHR) $(MFOBJ) $(MFCAP) $(NUMOBJ) $(QAG) $(WFUNC) $(RFUNC)
	${FC} -o gentest.x $(MAIN) $(MFSHR) $(MFOBJ) $(MFCAP) $(NUMOBJ) $(QAG) $(WFUNC) $(RFUNC)

gencaplib.so: $(MFSHR) $(MFOBJ) $(MFCAP) $(NUMOBJ) $(QAG) $(WFUNC) $(RFUNC)
	$(FC) $(FOPT) -shared -o $@ $(MFSHR) $(MFOBJ) $(MFCAP) $(NUMOBJ) $(QAG) $(WFUNC) $(RFUNC)

gputest.x: $(MAIN) $(MFSHR) $(MFOBJ) $(MFCAP) $(MFGPU) $(NUMOBJ) $(QAG) $(WFUNC) $(RFUNCDEV)
	${FC} -o gentest.x $(MAIN) $(MFSHR) $(MFOBJ) $(MFCAP) $(MFGPU) $(NUMOBJ) $(QAG) $(WFUNC) $(RFUNCDEV)


$(MAIN): %.o: %.f90
	$(FC) $(FOPT) -c  $<

$(MFSHR): %.o: %.f90
	$(FC) $(FOPT) -c  $<

$(MFOBJ): %.o: %.f90
	$(FC) $(FOPT) -c  $<

$(MFGPU): %.o: %.f90
	$(FC) -Mcuda $(FOPT) -c $<

$(MFCAP): %.o: %.f90
	$(FC) $(FOPT) -c  $<

$(NUMOBJ): %.o: $(NUMDIR)/%.f
	$(FC) $(FOPT) -c  $<

$(QAG): %.o: $(QAGDIR)/%.f
	$(FC) $(FOPT) -c  $<

$(WFUNC): %.o: $(WDIR)/%.f
	$(FC) $(FOPT) -c  $<

$(RFUNC): %.o: $(RDIR)/%.f
	$(FC) $(FOPT) -c  $<

$(RFUNCDEV): %.o: $(RDIRDEV)/%.f
	$(FC) -Mcuda $(FOPT) -c  $<

clean:
	rm -f *.o *.so *.mod gentest.x
