LIBS=`root-config --libs`
CFLAGS=`root-config --cflags` -c -g
LDFLAGS=`root-config --glibs`
CC=g++



#DalitzSmear: DalitzSmear.h DalitzSmear.cpp
#	$(CC) DalitzSmear.cpp $(LIBS) $(CFLAGS) -o smear

all: sm

sm: DalitzSmear.o mathUtility.o
	g++ DalitzSmear.o mathUtility.o $(LDFLAGS) $(LIBS) -o sm_dir

DalitzSmear.o: DalitzSmear.cpp DalitzSmear.h ../mathUtility/mathUtility.h
	g++ DalitzSmear.cpp $(CFLAGS)

mathUtility.o: ../mathUtility/mathUtility.cpp ../mathUtility/mathUtility.h ../vectorFactory/vectorFactory.h
	g++ ../mathUtility/mathUtility.cpp $(CFLAGS)

clean:
	rm *o sm_dir
