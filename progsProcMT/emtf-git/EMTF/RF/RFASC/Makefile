FC = f77
FFLAGS= -O
BIN_DIR =	../../test/bin

rfasc:	rfasc.f inpu_asc.f ../out_bin.f ../../include/datsz.inc nch.inc format.inc ../../include/four_byte.inc
	$(FC) $(FFLAGS) -o rfasc rfasc.f inpu_asc.f ../out_bin.f

clean:
	rm -f *.o

install:	
		\mv rfasc $(BIN_DIR)/.
