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

#run the simulation to create the event files: CM, unsmeared, and smeared .hepevt
echo "beginning simulation"
#this pipes output from sim into log for debugging
#./sim $Nevs $M $m_e $initP $chiCont > EventOutputs/log.txt
./sim $Nevs $M $m_e $initP $chiCont $scaleParameterNP

echo "simulation complete"
echo "beginning minimization"
#switch directories because MINUIT needs its own program 
cd Minimization

#run the numerical minimization and output the fit results to file, also pipes std:out into its own file so MINUIT doesn't print to screen N times
./min $Nevs $M $m_e $initP> ../EventOutputs/MinuitOutputDump.txt

echo "minimization complete"
echo "beginning plotting"
#switch directories over to the histogramming class and result interpreter for pull distributions
cd ..
cd plotTools/ResultInterpreter

#plot the pull distributions with this result interpreter program
./hist ../../EventOutputs/EventResults.txt $Nevs $initP ../../EventOutputs/ResultHistos.root

#switch over to the histogramming class
cd ..
cd Histogrammar

#run the program 3 times 1 for each mode plots should be the same, only input file stream is different and the boundaries for certain plots

./hist "../../EventOutputs/DalitzCMVectors.hepevt" $Nevs $US_CMplot $M $m_e $initP
./hist "../../EventOutputs/DalitzActualVectors.hepevt" $Nevs $US_LABplot $M $m_e $initP
./hist "../../EventOutputs/DalitzSmearVectors.hepevt" $Nevs $SM_LABplot $M $m_e $initP

echo "plotting complete"
