LIBS=`root-config --libs`
LDFLAGS=`root-config --glibs`
CFLAGS=`root-config --cflags` -c -g
CC=g++

#all : NumericalMinimization.cpp
#	$(CC) NumericalMinimization.cpp -o mins $(LIBS) $(CFLAGS)

all: thetastar

thetastar: MpiTheta.o MinHelper.o vectorFactory.o eventGenerator.o XYgenerator.o mathUtility.o
	g++ MpiTheta.o MinHelper.o vectorFactory.o eventGenerator.o XYgenerator.o mathUtility.o -o thetastar $(LIBS) 

MpiTheta.o: MpiTheta.cpp
	g++ MpiTheta.cpp $(CFLAGS)

MinHelper.o: ../../Minimization/MinHelper.cpp ../../Minimization/MinHelper.h ../../vectorFactory/vectorFactory.h #../Chisq/DalitzChiSq.h
	g++ ../../Minimization/MinHelper.cpp $(CFLAGS)

vectorFactory.o: ../../vectorFactory/vectorFactory.cpp ../../vectorFactory/vectorFactory.h ../../eventGenerator/eventGenerator.h
	g++ ../../vectorFactory/vectorFactory.cpp $(CFLAGS)

eventGenerator.o: ../../eventGenerator/eventGenerator.cpp ../../eventGenerator/eventGenerator.h ../../MC_Rejection/XYgenerator.h  #XYgenerator.o
	g++ ../../eventGenerator/eventGenerator.cpp $(CFLAGS)

XYgenerator.o: ../../MC_Rejection/XYgenerator.cpp ../../MC_Rejection/XYgenerator.h
	g++ ../../MC_Rejection/XYgenerator.cpp $(CFLAGS)

mathUtility.o: ../../mathUtility/mathUtility.cpp ../../mathUtility/mathUtility.h ../../vectorFactory/vectorFactory.h
	g++ ../../mathUtility/mathUtility.cpp $(CFLAGS)

clean:
	rm *o thetastar
