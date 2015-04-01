# ........................................................................
#f77 
	FF        = gfortran

#g95
#	FF        = g95
#.........................................................................
	FFLAGS    = -c -C
OBJINV = WALDIM.f inoutdata.f external.f

WALDIMALL : $(OBJINV)
	$(FF) $(OBJINV) -o WALDIM 


WALDIM: WALDIM.o inoutdata.o external.o
	$(FF) -o WALDIM WALDIM.o inoutdata.o external.o

WALDIM.o: WALDIM.f
	$(FF) $(FFLAGS) WALDIM.f
inoutdata.o: inoutdata.f
	$(FF) $(FFLAGS) inoutdata.f

external.o: external.f
	$(FF) $(FFLAGS) external.f
# .........................................................................
#
# f95
#OBJINV = WALDIM.f inoutdata.f external.f
#WALDIM: $(OBJINV)
#	f95 $(OBJINV) -o WALDIM 
# .........................................................................