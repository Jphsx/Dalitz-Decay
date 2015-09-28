LIBS=`root-config --libs`
CFLAGS=`root-config --cflags`
CC=g++



all: 
	$(CC) $(LIBS) $(CFLAGS) collectRMS.cpp -o rms
	$(CC) $(LIBS) $(CFLAGS) plotRMS.cpp -o plot

clean: 
	rm rms plot
