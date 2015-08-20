LIBS=`root-config --libs`
LDFLAGS=`root-config --glibs`
CFLAGS=`root-config --cflags` -c -g
CC=g++

#all : NumericalMinimization.cpp
#	$(CC) NumericalMinimization.cpp -o mins $(LIBS) $(CFLAGS)

all: min

min: MinHelper.o NumericalMinimization.o DalitzChiSq.o PhotonHelper.o vectorFactory.o eventGenerator.o XYgenerator.o mathUtility.o
	g++ MinHelper.o NumericalMinimization.o DalitzChiSq.o PhotonHelper.o vectorFactory.o eventGenerator.o XYgenerator.o mathUtility.o -o min $(LIBS) 

NumericalMinimization.o: NumericalMinimization.cpp MinHelper.h PhotonHelper.h ../mathUtility/mathUtility.h
	g++ NumericalMinimization.cpp $(CFLAGS)

MinHelper.o: MinHelper.cpp MinHelper.h ../vectorFactory/vectorFactory.h #../Chisq/DalitzChiSq.h
	g++ MinHelper.cpp $(CFLAGS)

DalitzChiSq.o: ../Chisq/DalitzChiSq.cpp ../Chisq/DalitzChiSq.h ../vectorFactory/vectorFactory.h ../mathUtility/mathUtility.h
	g++ ../Chisq/DalitzChiSq.cpp $(CFLAGS)

vectorFactory.o: ../vectorFactory/vectorFactory.cpp ../vectorFactory/vectorFactory.h ../eventGenerator/eventGenerator.h
	g++ ../vectorFactory/vectorFactory.cpp $(CFLAGS)

PhotonHelper.o: PhotonHelper.cpp PhotonHelper.h ../vectorFactory/vectorFactory.h ../mathUtility/mathUtility.h
	g++ PhotonHelper.cpp $(CFLAGS)

eventGenerator.o: ../eventGenerator/eventGenerator.cpp ../eventGenerator/eventGenerator.h ../MC_Rejection/XYgenerator.h  #XYgenerator.o
	g++ ../eventGenerator/eventGenerator.cpp $(CFLAGS)

XYgenerator.o: ../MC_Rejection/XYgenerator.cpp ../MC_Rejection/XYgenerator.h
	g++ ../MC_Rejection/XYgenerator.cpp $(CFLAGS)

mathUtility.o: ../mathUtility/mathUtility.cpp ../mathUtility/mathUtility.h ../vectorFactory/vectorFactory.h
	g++ ../mathUtility/mathUtility.cpp $(CFLAGS)

clean:
	rm *o min
