

LIBS=`root-config --libs`
CFLAGS=`root-config --cflags` -Wall -c -g
CC=g++

all: mc

mc: XYgenerator.o
	$(CC) XYgenerator.o $(LIBS) -o mc

XYgenerator.o: XYgenerator.cpp XYgenerator.h
	$(CC) XYgenerator.cpp $(LIBS) $(CFLAGS)

clean:
	rm *o mc
