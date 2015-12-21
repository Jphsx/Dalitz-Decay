#include "TFile.h"
#include <fstream>
#include <iomanip>
#include "TH1D.h"
#include <iostream>

int main(int argc, char* argv[]){

	double initP = atof(argv[1]);
	
	//let should add script to run all models, then individual rms can be extracted from the rms file

// Path to Histograms  	
TFile *f = new TFile("../../EventOutputs/ResultHistos_Minimalist.root");
TFile *f_2 = new TFile("../../EventOutputs/ResultHistos_Minimalist2.root");
	//grab the pi0 energy sum from the result histos file
	//histogram of the minimalist energy sum
	//TH1D * minRMS = (TH1D*)f->Get("hESumConstrained");
	//histogram of the elimination method and minimization fit Energy vaules
	//TH1D * elimRMS = (TH1D*)f->Get("hEfitSum");
	//measured value of pi0 energy
	//TH1D * measRMS = (TH1D*)f->Get("hESum");


	//pull the RMS from the histograms
	//double 	min = minRMS->GetRMS();
	//double elim = elimRMS->GetRMS();
	//double meas = measRMS->GetRMS();

//        double myrms = (TH1D*)f->Get("hEfitSum")->GetRMS();


	//get the minimimalist stuff
	TH1D * minRMS = (TH1D*)f->Get("hmass_min");
	TH1D * measRMS = (TH1D*)f->Get("hmass_m");
	TH1D * min2RMS = (TH1D*)f_2->Get("hmass_min2");

	double min = minRMS->GetRMS();
	double meas = measRMS->GetRMS();
	double min2 = min2RMS->GetRMS();
	//append the RMS to a text file
	
	std::ofstream f1;
	std::ofstream f2;
	f1.open("RMSmin.txt", std::ofstream::out | std::ofstream::app);
	f2.open("RMSmin2.txt", std::ofstream::out | std::ofstream::app);
	f1<<std::setprecision(9);
	f2<<std::setprecision(9);
	
	//f1<<meas<<" "<<min<<" "<<elim<<" "<<initP<<std::endl;
	f1<<meas<<" "<<min<<" "<<initP<<std::endl;
	f2<<meas<<" "<<min2<<" "<<initP<<std::endl;

	f1.close();
	f2.close();
	
	return 0;


}
