#!/bin/sh

rm EventOutputs/DalitzContour.root
rm EventOutputs/ResultHistos.root

M="0.13497"
m_e="0.000511"
initP="10.0"
Nevs="10000"
scaleParameterNP="1e-3" #leaves scale parameters to default in directional smearing
chiCont="0" #0:false 1:True
US_CMplot="0" #mode 0 plots CM event
US_LABplot="1" #mode 1 plots unsmeared Lab event
SM_LABplot="2" #mode 2 plots smeared Lab event
#if one of the histogramming programs is unnecessary, set variable to -1

#seed to control the RNG objects in the simulation
seed="0" #0 defaults to a random seed based on TRandom1

#run the simulation to create the event files: CM, unsmeared, and smeared .hepevt
echo "beginning simulation"
cd ../..
#this pipes output from sim into log for debugging
#./sim $Nevs $M $m_e $initP $chiCont > EventOutputs/log.txt
./sim $Nevs $M $m_e $initP $chiCont $scaleParameterNP $seed

echo "simulation complete"

cd scriptTools/Mpi_thetaStar

./thetastar $Nevs $M $m_e $initP
