#!/bin/sh


#maxSigma=.001
#minSigma=.000001
stepSize="25.0"


M="0.13497"
m_e="0.000511"
initP="10.0"
maxP="350.0"
Nevs="1000"
scaleParameterNP="0.0" #leaves scale parameters to default in directional smearing
chiCont="0" #0:false 1:True


rm RMS.txt #delete the old rms file so data doesnt get appended to it
#move up into the main directory
cd ../..

while [ $( echo "$initP < $maxP" | bc) -eq 1 ]
do

	#run the simulation to create the event files: CM, unsmeared, and 		smeared .hepevt
	echo "beginning simulation"
	#this pipes output from sim into log for debugging
	#./sim $Nevs $M $m_e $initP $chiCont > EventOutputs/log.txt
	./sim $Nevs $M $m_e $initP $chiCont $scaleParameterNP

	echo "simulation complete"
	echo "beginning minimization"
	#switch directories because MINUIT needs its own program 
	cd Minimization

	#run the numerical minimization and output the fit results to file, also 	pipes std:out into its own file so MINUIT doesn't print to screen N times
	./min $Nevs $M $m_e > ../EventOutputs/MinuitOutputDump.txt

	echo "minimization complete"
	echo "beginning plotting"
	#switch directories over to the histogramming class and result interpreter 	for pull distributions
	cd ..
	cd plotTools/ResultInterpreter

	#plot the pull distributions with this result interpreter program
	./hist ../../EventOutputs/EventResults.txt $Nevs

	cd ../../scriptTools/Epi_initP
	##ADD CODE--------
	# extract RMS from Epi0 plots (maybe write them to a file) then have a later script construct a histogram
	./rms $initP
	
	cd ../..

	#adjust the sigma value and rerun the simulation
	echo "Initial momentum set to $initP "
	initP=`echo "$initP + $stepSize" | bc`
	
done
cd scriptTools/Epi_initP
./plot
