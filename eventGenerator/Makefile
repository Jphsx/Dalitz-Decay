LIBS=`root-config --libs`
CFLAGS=`root-config --cflags` -Wall -c -g
CC=g++


all: gen

gen: eventGenerator.o XYgenerator.o Histogrammar.o
	$(CC) eventGenerator.o XYgenerator.o Histogrammar.o -o gen $(LIBS) 

eventGenerator.o: eventGenerator.cpp eventGenerator.h ../MC_Rejection/XYgenerator.h ../plotTools/Histogrammar/Histogrammar.h
	$(CC) eventGenerator.cpp $(LIBS) $(CFLAGS)

XYgenerator.o: ../MC_Rejection/XYgenerator.cpp ../MC_Rejection/XYgenerator.h
	$(CC) ../MC_Rejection/XYgenerator.cpp $(LIBS) $(CFLAGS)

Histogrammar.o: ../plotTools/Histogrammar/Histogrammar.h ../plotTools/Histogrammar/Histogrammar.cpp
	$(CC) ../plotTools/Histogrammar/Histogrammar.cpp $(LIBS) $(CFLAGS)

clean:
	rm *o gen
