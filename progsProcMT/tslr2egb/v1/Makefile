objects = main.o bytes1.o calendar.o declination.o startcount.o readtbl.o writesp.o writecal.o

ts2egb: $(objects)
        g++ -o ts2egb $(objects)

main.o: main.cpp bytes1.h count.h

bytes1.o: bytes1.cpp bytes1.h

calendar.o: calendar.c calendar.h

declination.o: declination.cpp declination.h

startcount.o: startcount.cpp count.h calendar.h bytes1.h

readtbl.o : readtbl.cpp readtbl.h count.h bytes1.h

writesp.o: writesp.cpp readtbl.h declination.h

writecal.o: writecal.cpp

clean:
       rm edit $(objects)
