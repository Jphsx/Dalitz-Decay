LIBS=`root-config --libs`
CFLAGS=`root-config --cflags` -Wall -g -c
CC=g++


all: vf

vf: vectorFactory.o eventGenerator.o XYgenerator.o
	$(CC) vectorFactory.o eventGenerator.o XYgenerator.o -o vf $(LIBS) 

vectorFactory.o: vectorFactory.cpp vectorFactory.h ../eventGenerator/eventGenerator.h #../MC_Rejection/XYgenerator.h
	$(CC) vectorFactory.cpp $(LIBS) $(CFLAGS)

eventGenerator.o: ../eventGenerator/eventGenerator.cpp ../eventGenerator/eventGenerator.h ../MC_Rejection/XYgenerator.h 
	$(CC) ../eventGenerator/eventGenerator.cpp $(LIBS) $(CFLAGS)

XYgenerator.o: ../MC_Rejection/XYgenerator.cpp ../MC_Rejection/XYgenerator.h
	$(CC) ../MC_Rejection/XYgenerator.cpp $(LIBS) $(CFLAGS)

clean:
	rm *o vf
