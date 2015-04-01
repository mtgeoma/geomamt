# Make tslr2egb

OBJS=		main.o bytes.o calendar.o declination.o readtbl.o startcount.o writecal.o writeegb.o writesp.o

clean:
		rm -f *.o *~ core

tslr2egb:		$(OBJS)
		g++ -o $@ $(OBJS)

main.o: main.cpp bytes.h count.h readtbl.h
	g++ -c main.cpp

bytes.o: bytes.cpp bytes.h
	g++ -c bytes.cpp

calendar.o: calendar.c calendar.h
	g++ -c calendar.c

declination.o: declination.cpp declination.h
	g++ -c declination.cpp

readtbl.o : readtbl.cpp readtbl.h
	g++ -c readtbl.cpp

startcount.o: startcount.cpp count.h
	g++ -c startcount.cpp

writecal.o: writecal.cpp readtbl.h
	g++ -c writecal.cpp

writeegb.o: writeegb.cpp bytes.h count.h
	g++ -c writeegb.cpp

writesp.o: writesp.cpp readtbl.h
	g++ -c writesp.cpp
