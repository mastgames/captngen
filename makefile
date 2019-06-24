FC=gfortran
#FOPT= -O3
FOPT= -O3 -fPIC
AUXDIR = ./aux
QAGDIR = ./aux/dqag
WDIR = ./Wfunctions
RDIR = ./Rfunctions

MAIN = main.o
MFOBJ = gencap.o
TRGOBJ = alphakappamod.o transgen.o fastevap.o
NUMOBJ =  dgamic.o d1mach.o
AUXOBJ = sgolay.o spline.o
QAG=  dsntdqagse.o dqelg.o dqk21.o dqpsrt.o dsntdqk21.o
WFUNC = WM.o WS2.o WS1.o WP2.o WMP2.o WP1.o WD.o WS1D.o
RFUNC = RM.o RS2.o RS1.o RP2.o RMP2.o RP1.o RD.o RS1D.o


gentest.x: $(MAIN) $(MFOBJ) $(TRGOBJ) $(NUMOBJ) $(AUXOBJ) $(QAG) $(WFUNC) $(RFUNC)
	${FC} -o gentest.x $(MAIN) $(MFOBJ) $(TRGOBJ) $(NUMOBJ) $(AUXOBJ) $(QAG) $(WFUNC) $(RFUNC)
#	rm $(MFOBJ) $(NUMOBJ) $(QAG)


gencaplib.so: $(MFOBJ) $(TRGOBJ) $(NUMOBJ) $(AUXOBJ) $(QAG) $(WFUNC) $(RFUNC)
	$(FC) $(FOPT) -shared -o $@ $(MFOBJ) $(TRGOBJ) $(NUMOBJ) $(AUXOBJ) $(QAG) $(WFUNC) $(RFUNC)


$(AUXOBJ): %.o : $(AUXDIR)/%.f90
	$(FC) $(FOPT) -c  $<

$(QAG): %.o : $(QAGDIR)/%.f
	$(FC) $(FOPT) -c  $<

$(MFOBJ): %.o: %.f90
	$(FC) $(FOPT) -c  $<

$(TRGOBJ): %.o: %.f90
	$(FC) $(FOPT) -c $<

$(MAIN): %.o: %.f90
	$(FC) $(FOPT) -c  $<

$(NUMOBJ): %.o: $(AUXDIR)/%.f
	$(FC) $(FOPT) -c  $<

$(WFUNC): %.o: $(WDIR)/%.f
	$(FC) $(FOPT) -c  $<

$(RFUNC): %.o: $(RDIR)/%.f
	$(FC) $(FOPT) -c  $<


clean:
	rm -f *.o *.so gentest.x
