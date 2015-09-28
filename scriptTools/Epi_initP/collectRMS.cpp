#include "TFile.h"
#include <fstream>
#include <iomanip>
#include "TH1D.h"
#include <iostream>

int main(int argc, char* argv[]){

	double initP = atof(argv[1]);

// Path to Histograms  	
TFile *f = new TFile("../../EventOutputs/ResultHistos.root");

	//grab the pi0 energy sum from the result histos file
	//histogram of the minimalist energy sum
	TH1D * minRMS = (TH1D*)f->Get("hESumConstrained");
	//histogram of the elimination method and minimization fit Energy vaules
	TH1D * elimRMS = (TH1D*)f->Get("hEfitSum");
	//measured value of pi0 energy
	TH1D * measRMS = (TH1D*)f->Get("hESum");

	//pull the RMS from the histograms
	double 	min = minRMS->GetRMS();
	double elim = elimRMS->GetRMS();
	double meas = measRMS->GetRMS();

//        double myrms = (TH1D*)f->Get("hEfitSum")->GetRMS();

	//append the RMS to a text file
	
	std::ofstream f1;
	f1.open("RMS.txt", std::ofstream::out | std::ofstream::app);
	
	f1<<std::setprecision(9);
	
	f1<<meas<<" "<<min<<" "<<elim<<" "<<initP<<std::endl;
	
	f1.close();
	
	return 0;


}
