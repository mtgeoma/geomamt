#	Makefile for DNFF
#	edit the BIN_DIR directory to install dnff in your binary directory
#	use make all to compile, link, install dnff and remove all .o files
#	for debugging exchange the -O2 and -g options in FFLAGS


FC = gfortran
FFLAGS= -O2 -Wall -mcmodel=medium
# -g
BIN_DIR= ../test/bin


TARGETS=	dnff install clean

OBJS=		afcor.o autocor.o badrec.o bdrcsu.o ft_subs.o chdec.o \
		cldrft.o cmove.o dcfilt.o dcimte.o decset.o demean.o \
		dnff.o dtrnd.o fcorsu.o filtcor.o freqout.o frstdif.o \
		getsp.o ltslv.o mkset.o ptdist.o pterst.o sort.o \
		resptbl.o inpu_bin.o out_pack.o mk_offst.o phs_shft.o \
		inpu_asc.o metronix.o fitpack.o ssort.o str_utils.o

all:		$(TARGETS)

clean:		
		rm -f *.o *~ core

install:
		\mv dnff $(BIN_DIR)/.

dnff:		$(OBJS)
		$(FC) $(FFLAGS) -o $@ $(OBJS)
