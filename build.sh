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

cd ../..
cd Minimalist
make clean
make

cd ..
cd Minimalist2
make clean
make

cd ..
cd scriptTools

cd Epi_initP
make clean
make
cd ..

cd Epi_SP
make clean
make
cd ..

cd Mpi_initP
make clean
make
cd ..

cd Sigma_Psi123
make clean
make
cd ../.. 


