LIBS=`root-config --libs`
LDFLAGS=`root-config --glibs`
CFLAGS=`root-config --cflags` -c -g
CC=g++

#all : DalitzChiSq2.cpp
#	$(CC) DalitzChiSq2.cpp -o chisq $(LIBS) $(CFLAGS)

all: sim

sim: DalitzSimulation.o DalitzSmear.o DalitzChiSq.o vectorFactory.o eventGenerator.o XYgenerator.o mathUtility.o
	g++ DalitzSimulation.o DalitzSmear.o DalitzChiSq.o vectorFactory.o eventGenerator.o XYgenerator.o mathUtility.o $(LDFLAGS) $(LIBS) -o sim

DalitzSimulation.o: Simulation/DalitzSimulation.cpp Smearing/DalitzSmear.h Chisq/DalitzChiSq.h mathUtility/mathUtility.h
	g++ Simulation/DalitzSimulation.cpp $(CFLAGS)

DalitzSmear.o: Smearing/DalitzSmear.cpp Smearing/DalitzSmear.h mathUtility/mathUtility.h
	g++ Smearing/DalitzSmear.cpp $(CFLAGS)

DalitzChiSq.o: Chisq/DalitzChiSq.cpp Chisq/DalitzChiSq.h vectorFactory/vectorFactory.h mathUtility/mathUtility.h
	g++ Chisq/DalitzChiSq.cpp $(CFLAGS)

vectorFactory.o: vectorFactory/vectorFactory.cpp vectorFactory/vectorFactory.h eventGenerator/eventGenerator.h
	g++ vectorFactory/vectorFactory.cpp $(CFLAGS)

eventGenerator.o: eventGenerator/eventGenerator.cpp eventGenerator/eventGenerator.h MC_Rejection/XYgenerator.h  #XYgenerator.o
	g++ eventGenerator/eventGenerator.cpp $(CFLAGS)

XYgenerator.o: MC_Rejection/XYgenerator.cpp MC_Rejection/XYgenerator.h
	g++ MC_Rejection/XYgenerator.cpp $(CFLAGS)

mathUtility.o: mathUtility/mathUtility.cpp mathUtility/mathUtility.h vectorFactory/vectorFactory.h
	g++ mathUtility/mathUtility.cpp $(CFLAGS)

clean:
	rm *o sim
