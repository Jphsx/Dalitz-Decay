#!/bin/sh


#maxSigma=.001
#minSigma=.000001
stepSize="1"

M="0.13497"
m_e="0.000511"
initP="1"
maxP="75"

#testing for minimalist 2
#stepSize="2"
#initP="10"
#maxP="11"

##############

Nevs="1000"
scaleParameterNP="1e-4" #leaves scale parameters to default in directional smearing
chiCont="0" #0:false 1:True


rm RMSmin.txt #delete the old rms file so data doesnt get appended to it
rm RMSmin2.txt
#move up into the main directory
cd ../..
cd EventOutputs
rm ResultHistos_Minimalist.root
rm ResultHistos_Minimalist2.root
cd ..

while [ $( echo "$initP < $maxP" | bc) -eq 1 ]
do

	#run the simulation to create the event files: CM, unsmeared, and 		smeared .hepevt
	echo "beginning simulation"
	#this pipes output from sim into log for debugging
	#./sim $Nevs $M $m_e $initP $chiCont > EventOutputs/log.txt
	./sim $Nevs $M $m_e $initP $chiCont $scaleParameterNP

	echo "simulation complete"
	echo "beginning minimization"
	
	cd Minimalist2
	./min2 $Nevs $M $m_e $initP 
	echo "end min2"

	cd ..
	#switch directories because MINUIT needs its own program 
	#cd Minimization #remember minimization doesnt work because minuit wont converge, photon is too smeared
	cd Minimalist
	#run the numerical minimization and output the fit results to file, also 	pipes std:out into its own file so MINUIT doesn't print to screen N times
	./min $Nevs $M $m_e $initP

	echo "minimization complete"
	echo "beginning plotting"
	#switch directories over to the histogramming class and result interpreter 	for pull distributions
	cd ..
	cd plotTools/ResultInterpreter

	#plot the pull distributions with this result interpreter program
	./hist ../../EventOutputs/EventResults_Minimalist.txt $Nevs $initP ../../EventOutputs/ResultHistos_Minimalist.root
	./hist ../../EventOutputs/EventResults_Minimalist2.txt $Nevs $initP ../../EventOutputs/ResultHistos_Minimalist2.root	
	echo "end plotting"
	echo "begin scriptTool"
	
	cd ../../scriptTools/Epi_initP
	##ADD CODE-------- (to rms and plot for minimalist 2)
	# extract RMS from Epi0 plots (maybe write them to a file) then have a later script construct a histogram
	./rms $initP

	
	cd ../..

	#adjust the sigma value and rerun the simulation
	echo "Initial momentum set to $initP "
	initP=`echo "$initP + $stepSize" | bc`
	
done
cd scriptTools/Epi_initP
./plot
