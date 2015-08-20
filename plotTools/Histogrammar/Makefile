LIBS=`root-config --libs`
CFLAGS=`root-config --cflags` -c -g -Wall
CC=g++


all: hist

hist: Histogrammar.o
	$(CC) Histogrammar.o -o hist $(LIBS)

Histogrammar.o: Histogrammar.h Histogrammar.cpp 
	$(CC) $(CFLAGS) Histogrammar.cpp

clean:
	rm *o hist
