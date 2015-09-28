#!/bin/sh

find . -type f -name '*.o' -delete
find . -type f -name '*~' -delete
make clean

cd eventGenerator
make clean
cd ..

cd mathUtility
make clean
cd ..

cd Minimization
make clean
cd ..

cd plotTools/Histogrammar
make clean
cd ..

cd ResultInterpreter
make clean
cd ../..

cd Simulation
make clean
cd ..

cd vectorFactory
make clean
cd ..

cd Chisq
make clean
cd ..

cd MC_Rejection
make clean
cd ..

cd Smearing
make clean
cd ..

cd scriptTools/Epi_SP
make clean
cd ..