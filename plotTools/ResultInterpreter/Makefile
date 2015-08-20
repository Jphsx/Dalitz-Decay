LIBS=`root-config --libs`
LDFLAGS=`root-config --glibs`
CFLAGS=`root-config --cflags` -Wall -c -g
CC=g++

all: hist

hist: ResultInterpreter.o
	$(CC) ResultInterpreter.o -o hist $(LIBS)

ResultInterpreter.o: ResultInterpreter.h ResultInterpreter.cpp
	$(CC) ResultInterpreter.cpp $(CFLAGS)

clean:
	rm *o hist
