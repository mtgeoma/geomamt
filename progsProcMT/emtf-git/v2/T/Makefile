#	Makefile for TRANMT
#
FC = gfortran
FFLAGS= -O2 -Wall
# -g
BIN_DIR= ../test/bin

TARGETS=	tranmt install clean

OBJS=		bset.o chdec.o cmove.o edfwts.o matmult.o mkfdir.o \
		mkhat.o mkrec.o mkrhat.o modobs.o movec.o mult.o \
		rbstrn.o rdcndwt.o rfhead.o mcohc.o rxspclev.o \
		savtf.o setup.o sindex.o sort.o stack.o trlsrr.o \
		unstack.o wt.o sdmsort.o cohsrt.o cohband.o \
		rtreref.o tranmt.o fop.o pwrvf.o readfvg.o wrt_z.o \
		reorder.o sdmwrite.o corchng.o rotsdm.o

clean:
	rm -f *.o core

install:
		\mv tranmt $(BIN_DIR)/.

tranmt:		$(OBJS)
		$(FC) $(FFLAGS) -o $@ $(OBJS)
