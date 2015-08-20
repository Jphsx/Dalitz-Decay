LIBS=`root-config --libs`
CFLAGS=`root-config --cflags` -Wall -c 
CC=g++

#all : DalitzSimulation.cpp ../Smearing/DalitzSmear.h ../vectorFactory/vectorFactory.h ../Chisq/DalitzChiSq.h 
#	$(CC) DalitzSimulation.cpp -o sim $(LIBS) $(CFLAGS)

all: sim

sim: DalitzSimulation.o DalitzSmear.o DalitzChiSq.o vectorFactory.o
	g++ DalitzSimulation.o DalitzSmear.o DalitzChiSq.o vectorFactory.o -o sim $(LIBS) #$(CFLAGS)

DalitzSimulation.o: DalitzSimulation.cpp
	g++ DalitzSimulation.cpp $(LIBS) $(CFLAGS)

DalitzSmear.o: ../Smearing/DalitzSmear.cpp ../Smearing/DalitzSmear.h
	g++ ../Smearing/DalitzSmear.cpp $(LIBS) $(CFLAGS)

DalitzChiSq.o: ../Chisq/DalitzChiSq.cpp ../Chisq/DalitzChiSq.h
	g++ ../Chisq/DalitzChiSq.cpp $(LIBS) $(CFLAGS)

vectorFactory.o: ../vectorFactory/vectorFactory.cpp ../vectorFactory/vectorFactory.h 
	g++ ../vectorFactory/vectorFactory.cpp $(LIBS) $(CFLAGS)

clean:
	rm *o sim
