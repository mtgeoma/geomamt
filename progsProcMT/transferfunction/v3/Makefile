# Make datetext

OBJS=		main.o transferfunction.o commandline.o inputdata.o emprogformat.o ediformat.o pt1format.o egbertformat.o jonesformat.o

clean:
		rm -f *.o *~ core

transferfunction:	$(OBJS)
		g++ -o $@ $(OBJS)

main.o: main.cpp transferfunction.h commandline.h inputdata.h
	g++ -c main.cpp

transferfunction.o: transferfunction.cpp transferfunction.h
	g++ -c transferfunction.cpp

commandline.o: commandline.cpp commandline.h
	g++ -c commandline.cpp

inputdata.o: inputdata.cpp inputdata.h
	g++ -c inputdata.cpp

emprogformat.o: emprogformat.cpp inputdata.h
	g++ -c emprogformat.cpp

ediformat.o: ediformat.cpp inputdata.h
	g++ -c ediformat.cpp

pt1format.o: pt1format.cpp inputdata.h
	g++ -c pt1format.cpp

egbertformat.o: egbertformat.cpp inputdata.h
	g++ -c egbertformat.cpp

jonesformat.o: jonesformat.cpp inputdata.h
	g++ -c jonesformat.cpp
