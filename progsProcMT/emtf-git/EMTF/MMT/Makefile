##  Multiple station processing program

#  change BIN_DIR to desired bin directory 
BIN_DIR= ../test/bin
FC = gfortran

#  options for debug compilation
FCFLAGS =

#  compiler options for Sun Solaris
FFLAGS = -O2

#  following are options for compiling on IBM RS 6000/SP2 cluster
#FFLAGS = -O -qrealsize=4 -qfixed -qmaxmem=-1 \
#         -L/usr/local/lib -bmaxdata:512000000 \
#         -bmaxstack:512000000 -qSPILL=3270
# end of IBM options

#  Libraries to link to (assuming liblapack.a and liblas.a are in
#				/usr/local/lib  )
LIBS = -llapack -lblas

### ALTERNATIVE :: local, very partial lapack library
###    that you can make with source code include with EMTF/MMT
###    (only QR, general linear solver, BLAS)
###FFLAGS = -O -L../lapack

TARGETS=	multmtrn install clean

OBJS=		anfld.o mk_fd_rec.o mmt_mtrx.o mmt_subs.o multmtrn.o \
		prteig.o rbstk2.o readfvg.o rsp.o setup.o timerr.o \
		ln_rbst.o sindex.o extra.o n_rbst.o grad.o rbstreg.o \
		var_adj.o rr_rbst.o pc_out.o wrtx.o \
		filtpc.o canonic.o cn_wt.o sep_s_n.o wrt_z.o \
		refproj.o mk_list.o corchng.o minpwr.o

all:		$(TARGETS)

clean:
		rm -f *.o *~ core


install:
		\mv multmtrn $(BIN_DIR)/.

multmtrn:	$(OBJS)
		$(FC) $(FFLAGS) -o $@ $(OBJS) $(LIBS)
