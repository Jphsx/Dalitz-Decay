LIBS=`root-config --libs`
LDFLAGS=`root-config --glibs`
CFLAGS=`root-config --cflags` -c -g
CC=g++


all: xy

xy: XY_surf.o XYgenerator.o
	g++ XY_surf.o XYgenerator.o -o xy $(LIBS) 

XYgenerator.o: ../../MC_Rejection/XYgenerator.cpp ../../MC_Rejection/XYgenerator.h
	g++ ../../MC_Rejection/XYgenerator.cpp $(CFLAGS)

XY_surf.o: XY_surf.cpp
	g++ XY_surf.cpp $(CFLAGS)

clean:
	rm *o xy
 
