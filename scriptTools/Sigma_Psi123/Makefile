LIBS=`root-config --libs`
LDFLAGS=`root-config --glibs`
CFLAGS=`root-config --cflags` -c -g
CC=g++

#all : NumericalMinimization.cpp
#	$(CC) NumericalMinimization.cpp -o mins $(LIBS) $(CFLAGS)

all: angle

angle: Psi123RMS.o vectorFactory.o mathUtility.o eventGenerator.o XYgenerator.o DalitzSmear.o
	g++ Psi123RMS.o vectorFactory.o mathUtility.o eventGenerator.o XYgenerator.o DalitzSmear.o -o angle $(LIBS) 

Psi123RMS.o: Psi123RMS.cpp ../../Smearing/DalitzSmear.h ../../mathUtility/mathUtility.h
	g++ Psi123RMS.cpp $(CFLAGS)

DalitzSmear.o: ../../Smearing/DalitzSmear.h ../../Smearing/DalitzSmear.cpp
	g++ ../../Smearing/DalitzSmear.cpp $(CFLAGS)

eventGenerator.o: ../../eventGenerator/eventGenerator.cpp ../../eventGenerator/eventGenerator.h ../../MC_Rejection/XYgenerator.h  #XYgenerator.o
	g++ ../../eventGenerator/eventGenerator.cpp $(CFLAGS)

XYgenerator.o: ../../MC_Rejection/XYgenerator.cpp ../../MC_Rejection/XYgenerator.h
	g++ ../../MC_Rejection/XYgenerator.cpp $(CFLAGS)


vectorFactory.o: ../../vectorFactory/vectorFactory.cpp ../../vectorFactory/vectorFactory.h 
	g++ ../../vectorFactory/vectorFactory.cpp $(CFLAGS)



mathUtility.o: ../../mathUtility/mathUtility.cpp ../../mathUtility/mathUtility.h ../../vectorFactory/vectorFactory.h
	g++ ../../mathUtility/mathUtility.cpp $(CFLAGS)



clean:
	rm *o angle
