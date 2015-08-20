LIBS=`root-config --libs`
LDFLAGS=`root-config --glibs`
CFLAGS=`root-config --cflags` -Wall -c -g 
CC=g++

#all : DalitzSimulation.cpp ../Smearing/DalitzSmear.h ../vectorFactory/vectorFactory.h ../Chisq/DalitzChiSq.h 
#	$(CC) DalitzSimulation.cpp -o sim $(LIBS) $(CFLAGS)

all: sim

sim: DalitzSimulation.o DalitzSmear.o DalitzChiSq.o vectorFactory.o
	g++ DalitzSimulation.o DalitzSmear.o DalitzChiSq.o vectorFactory.o $(LDFLAGS) $(LIBS)  -o sim

DalitzSimulation.o: DalitzSimulation.cpp DalitzSmear.h DalitzChiSq.h
	g++ DalitzSimulation.cpp $(LIBS) $(CFLAGS)

DalitzSmear.o: DalitzSmear.cpp DalitzSmear.h
	g++ DalitzSmear.cpp $(LIBS) $(CFLAGS)

DalitzChiSq.o: DalitzChiSq.cpp DalitzChiSq.h vectorFactory.h
	g++ DalitzChiSq.cpp $(LIBS) $(CFLAGS)

vectorFactory.o: vectorFactory.cpp vectorFactory.h
	g++ vectorFactory.cpp $(LIBS) $(CFLAGS)

#all: DalitzSimulation.cpp vectorFactory.h DalitzChiSq.h DalitzSmear.h

#	g++ DalitzSimulation.cpp $(LIBS) $(CFLAGS) 
clean:
	rm *o sim
