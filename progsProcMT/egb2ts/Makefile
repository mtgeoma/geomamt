# Make egb2ts

OBJSD=		main.o getconv.o writetsd.o
OBJSS=		main.o getconv.o writetss.o
OBJSISO=	main.o getconv.o writetsISO.o

clean:
		rm -f *.o *~ core

egb2tsd:		$(OBJSD)
		g++ -o $@ $(OBJSD)

egb2tss:		$(OBJSS)
		g++ -o $@ $(OBJSS)

egb2tsISO:		$(OBJSISO)
		g++ -o $@ $(OBJSISO) -lboost_date_time

main.o: main.cpp
	g++ -c main.cpp

getconv.o: getconv.cpp
	g++ -c getconv.cpp

writetsd.o: writetsd.cpp calendar.h
	g++ -c writetsd.cpp

writetss.o: writetss.cpp calendar.h
	g++ -c writetss.cpp

writetsISO.o: writetsISO.cpp
	g++ -c writetsISO.cpp
