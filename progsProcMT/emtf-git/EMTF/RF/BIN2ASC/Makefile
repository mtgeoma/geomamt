FC = f77
FFLAGS= -O
BIN_DIR =	../../test/bin
TARGETS=	bin2asc clean install

bin2asc:	bin2asc.f out_asc.f ../../include/datsz.inc input.f
		$(FC) $(FFLAGS) -o bin2asc bin2asc.f out_asc.f input.f

clean:
		rm -f *.o

install:	
		\mv bin2asc $(BIN_DIR)/.
