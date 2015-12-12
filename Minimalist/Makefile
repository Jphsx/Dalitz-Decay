LIBS=`root-config --libs`
LDFLAGS=`root-config --glibs`
CFLAGS=`root-config --cflags` -c -g
CC=g++

#all : NumericalMinimization.cpp
#	$(CC) NumericalMinimization.cpp -o mins $(LIBS) $(CFLAGS)

all: min

min: MinHelper.o minimalist.o DalitzChiSq.o vectorFactory.o eventGenerator.o XYgenerator.o mathUtility.o PhotonHelper.o
	g++ MinHelper.o minimalist.o DalitzChiSq.o vectorFactory.o eventGenerator.o XYgenerator.o mathUtility.o PhotonHelper.o -o min $(LIBS) 

minimalist.o: minimalist.cpp ../Minimization/MinHelper.h ../mathUtility/mathUtility.h
	g++ minimalist.cpp $(CFLAGS)

MinHelper.o: ../Minimization/MinHelper.cpp ../Minimization/MinHelper.h ../vectorFactory/vectorFactory.h #../Chisq/DalitzChiSq.h
	g++ ../Minimization/MinHelper.cpp $(CFLAGS)

DalitzChiSq.o: ../Chisq/DalitzChiSq.cpp ../Chisq/DalitzChiSq.h ../vectorFactory/vectorFactory.h ../mathUtility/mathUtility.h
	g++ ../Chisq/DalitzChiSq.cpp $(CFLAGS)

vectorFactory.o: ../vectorFactory/vectorFactory.cpp ../vectorFactory/vectorFactory.h ../eventGenerator/eventGenerator.h
	g++ ../vectorFactory/vectorFactory.cpp $(CFLAGS)

eventGenerator.o: ../eventGenerator/eventGenerator.cpp ../eventGenerator/eventGenerator.h ../MC_Rejection/XYgenerator.h  #XYgenerator.o
	g++ ../eventGenerator/eventGenerator.cpp $(CFLAGS)

XYgenerator.o: ../MC_Rejection/XYgenerator.cpp ../MC_Rejection/XYgenerator.h
	g++ ../MC_Rejection/XYgenerator.cpp $(CFLAGS)

mathUtility.o: ../mathUtility/mathUtility.cpp ../mathUtility/mathUtility.h ../vectorFactory/vectorFactory.h
	g++ ../mathUtility/mathUtility.cpp $(CFLAGS)

PhotonHelper.o: ../Minimization/PhotonHelper.cpp ../Minimization/PhotonHelper.h ../vectorFactory/vectorFactory.h #../Chisq/DalitzChiSq.h
	g++ ../Minimization/PhotonHelper.cpp $(CFLAGS)

clean:
	rm *o min
