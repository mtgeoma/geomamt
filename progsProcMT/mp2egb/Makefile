# Make mp2egb

OBJS=		main.o calendar.o declination.o meu_geomag70.o startcount.o readhead.o writesp.o hardware.o writeegb.o

clean:
		rm -f *.o *~ core

mp2egb:		$(OBJS)
		g++ -o $@ $(OBJS)

main.o: main.cpp count.h hardware.h
	g++ -c main.cpp

calendar.o: calendar.c calendar.h
	g++ -c calendar.c

declination.o: declination.cpp declination.h
	g++ -c declination.cpp

meu_geomag60.o: meu_geomag60.c meu_geomag60.h
	g++ -c meu_geomag60.c

startcount.o: startcount.cpp count.h
	g++ -c startcount.cpp

readhead.o : readhead.cpp readhead.h
	g++ -c readhead.cpp

writesp.o: writesp.cpp readhead.h hardware.h
	g++ -c writesp.cpp

hardware.o: hardware.cpp
	g++ -c hardware.cpp

writeegb.o: writeegb.cpp count.h
	g++ -c writeegb.cpp
