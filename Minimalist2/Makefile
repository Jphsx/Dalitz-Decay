LIBS=`root-config --libs`
LDFLAGS=`root-config --glibs`
CFLAGS=`root-config --cflags` -c -g
CC=g++

#all : NumericalMinimization.cpp
#	$(CC) NumericalMinimization.cpp -o mins $(LIBS) $(CFLAGS)

#make for helper testing
#helper: min2utility.o mathUtility.o
#	g++ min2utility.o mathUtility.o -o minhelper $(LIBS)

all: min2

min2: MinHelper.o minimalist2.o DalitzChiSq.o vectorFactory.o eventGenerator.o XYgenerator.o mathUtility.o PhotonHelper.o min2utility.o
	g++ MinHelper.o minimalist2.o DalitzChiSq.o vectorFactory.o eventGenerator.o XYgenerator.o mathUtility.o PhotonHelper.o min2utility.o -o min2 $(LIBS) 

minimalist2.o: minimalist2.cpp ../Minimization/MinHelper.h ../mathUtility/mathUtility.h min2utility.h
	g++ minimalist2.cpp $(CFLAGS)

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

min2utility.o: min2utility.cpp min2utility.h ../mathUtility/mathUtility.h
	g++ min2utility.cpp $(CFLAGS)

clean:
	rm *o min2 minhelper
