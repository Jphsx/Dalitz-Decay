#!/bin/sh

make clean
make
cd Minimization
make clean
make 
cd ..
cd plotTools/ResultInterpreter
make clean
make

cd ..
cd Histogrammar
make clean
make
