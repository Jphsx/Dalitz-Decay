#ifndef MINIMIZATION_
#define MINIMIZATION_
#include "MinHelper.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <iostream>
#include "TMath.h"
#include <cmath>
#include <stdio.h>
#include <cstdio>
#include "TLorentzVector.h"
#include "PhotonHelper.h"
#include <fstream>
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "../mathUtility/mathUtility.h"
#include <iomanip>


using namespace std;

/* Collection of all the variables associated with calculating fit values including MINUIT settings and output variables */
double x1m;
double variance1;
double theta1;

double x2m;
double variance2;
double theta2;

double x3m;
double variance3;
double M;
double m_e;
double z12;
double z23;
double z13;

double centerX1;
double centerX2;
double offset1;
double offset2;

double x1Fit;
double x2Fit;
double x1FitVariance;
double x2FitVariance;
double x3FitVariance;
double x12covariance;
double minChisq;

//numerical convergence 
double x3FitOld=-1;
double x3Fit=-1;
//double oldvariance;
double x1FitOld;
double x2FitOld;

double x1g,x2g,x3g;
double dpsi1,dpsi2,dpsi3;

     

//Algorthim that is to be minimized
double DalitzElimination(const double *vars){
		
	double term1 = (pow(vars[0]-x1m,2) )/ variance1;
	double term2 = (pow(vars[1]-x2m,2) )/ variance2;	

	//now for e3
	double topterm =   M*M - 2*m_e*m_e - (z12/(sin(theta1) * sin(theta2) *vars[0]*vars[1]));
	double bottomterm = ( (z13/(sin(theta1) * vars[0])) + (z23/(sin(theta2) * vars[1])) );
	double x3=topterm/bottomterm;

	double term3 = pow(x3-x3m,2)/variance3;
	

	return term1+term2+term3;

}
//function that uses MINUIT to perform the x^2 minimization
int NumericalMinimization(const char * minName = "Minuit2",
                          const char *algoName = "")
{

   ROOT::Math::Minimizer* min =
      ROOT::Math::Factory::CreateMinimizer(minName, algoName);

   // set tolerance , etc...
   min->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
   min->SetMaxIterations(10000000);  // for GSL
   min->SetTolerance(0.00000001);
   min->SetPrintLevel(1);

   // create funciton wrapper for minmizer
   // a IMultiGenFunction type
   ROOT::Math::Functor f(&DalitzElimination,2);
   double step[2] = {0.000001,0.000001};
   // starting point

   double variable[2] = { centerX1-offset1,centerX2-offset2};

   min->SetFunction(f);

   // Set the free variables to be minimized!
   min->SetVariable(0,"x",variable[0], step[0]);
   min->SetVariable(1,"y",variable[1], step[1]);

   // do the minimization
   min->Minimize();

   const double *xs = min->X();
   std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): "
             << min->MinValue()  << std::endl;

   if ( min->MinValue()  < 10  && f(xs) < 10)
      std::cout << "Minimizer " << minName << " - " << algoName
                << "   converged to the right minimum" << std::endl;
   else {
      std::cout << "Minimizer " << minName << " - " << algoName
                << "   Large X^2 value !!! " <<min->MinValue()<< std::endl;
     // Error("NumericalMinimization","fail to converge");
	}

	minChisq = min->MinValue();
	x1Fit=xs[0];
	x2Fit=xs[1];
	x1FitVariance=min->CovMatrix(0,0);
	x2FitVariance=min->CovMatrix(1,1);
	x12covariance=min->CovMatrix(0,1);
	//cout<<"COV"<< min->CovMatrix(0,1)<<endl;
	//cout<< min->CovMatrix(1,0)<<endl;
	std::cout<<"Variance estimate on x1_m "<<variance1<<std::endl;
	std::cout<<"Variance estimate on x2_m "<<variance2<<std::endl;
	std::cout<<"Variance estimate on x3_m "<<variance3<<std::endl;
	std::cout<<"Error Matrix on Fitted Parameters" << sqrt(x1FitVariance) << " " 
                 << sqrt(x2FitVariance) << " " << x12covariance/sqrt(x1FitVariance*x2FitVariance) << std::endl << std::endl;

   return 0;
}
//debugging method
void PrintVector(TLorentzVector v){
	cout<<v.Px()<<" "<<v.Py()<<" "<<v.Pz()<<" "<<v.P()<<" "<<v.E()<<" "<<v.M()<<endl;
}
int main(int argc, char* argv[]){
	//sets the number of minimzations to run i.e. N events
	int N=atoi(argv[1]);
	//invariant mass of the pion
	 M=atof(argv[2]);
	//electron mass
 	 m_e=atof(argv[3]);

	//initial momentum
	double initP=atof(argv[4]);
	
	
	TFile *froot = new TFile("../EventOutputs/ResultHistos.root", "UPDATE");
	/*TH1D* hESum = new TH1D("hESum", "Measured Energy Sum E1+E2+E3;pi0 Energy; Event per bin",100,9,11);
	TH1D* hEfitSum = new TH1D("hEfitSum", "Energy Sum from Minimization E1+E2+E3;pi0 Energy;Event per bin",100,9.92,10.07);
	TH1D* hESumConstrained = new TH1D("hESumConstrained", "Minimimalist E1+E2+E3C;pi0 Energy; event per bin",100,9.92,10.07);*/

	TH1D* hESum = new TH1D("hESum", "Measured Energy Sum E1+E2+E3;pi0 Energy; Event per bin",100,initP-0.5*initP,initP+0.5*initP);
	TH1D* hEfitSum = new TH1D("hEfitSum", "Energy Sum from Minimization E1+E2+E3;pi0 Energy;Event per bin",100,initP-0.5*initP,initP+0.5*initP);
	TH1D* hESumConstrained = new TH1D("hESumConstrained", "Minimimalist E1+E2+E3C;pi0 Energy; event per bin",100,initP-0.5*initP,initP+0.5*initP);
	
	TH1D* hEpiPull = new TH1D("hEpiPull", "pull distribution of pi0 energies; Epi_fit - Epi_m / sqrt(vfit-vm); event per bin",100,-0.5*initP,0.5*initP);
	double EpiPull;

	TH1D* hdpsix1 = new TH1D("hdpsix1","#delta#psi between generator and detector x1;Events per bin",100,-0.0,0.1);
	TH1D* hdpsix2 = new TH1D("hdpsix2","#delta#psi between generator and detector x2;Events per bin",100,-0.0,0.1);
	TH1D* hdpsix3 = new TH1D("hdpsix3","#delta#pis between generator and detector x3;Events per bin",100,-0.0,0.1);
	
	//each time the program is run the output which holds the fit values and uncertainties is wiped
	ofstream cleaner("../EventOutputs/EventResults.txt");
	cleaner<<"";
	cleaner.close();

		//reads and holds the particle structs to be set globally later
	MinHelper* h = new MinHelper("../EventOutputs/DalitzSmearVectors.hepevt","../EventOutputs/DalitzActualVectors.hepevt");

	//TTREE initiliaztion///////////
	TTree* tree = new TTree( "DalitzDecay" , "DalitzDecay");

	
	int evnum;
	//TLorentzVector epg_4v, emg_4v, gmg_4v, epm_4v, emm_4v, gmm_4v;
	TLorentzVector* genVects = new TLorentzVector[4]; //={epg_4v, emg_4v, gmg_4v};
	TLorentzVector* mesVects = new TLorentzVector[4]; //={epm_4v, emm_4v, gmm_4v};
	//TLorentzVector epg_4v= h->evtP_actual[1].v;
/*
	for(int i=1; i<=3; i++){
		genVects[i]=h->evtP_actual[i].v;
		mesVects[i]=h->evtP[i].v;
	}
*/	
        
        tree->Branch("evnum", &evnum);

        tree->Branch("x1m", &x1m);
        tree->Branch("x2m", &x2m);
        tree->Branch("x3m", &x3m);
        tree->Branch("x1f", &x1Fit);
        tree->Branch("x2f", &x2Fit);
        tree->Branch("x3f", &x3Fit);
        tree->Branch("x1g", &x1g);
        tree->Branch("x2g", &x2g);
        tree->Branch("x3g", &x3g);
	tree->Branch("x1mV", &variance1);
	tree->Branch("x2mV", &variance2);
	tree->Branch("x3mV", &variance3);
	tree->Branch("x1fV", &x1FitVariance);
	tree->Branch("x2fV", &x2FitVariance);
	tree->Branch("x3fV", &x3FitVariance);
	tree->Branch("p1g.", &genVects[1]);
	tree->Branch("p2g.", &genVects[2]);
	tree->Branch("p3g.", &genVects[3]);
	tree->Branch("p1m.", &mesVects[1]);
	tree->Branch("p2m.", &mesVects[2]);
	tree->Branch("p3m.", &mesVects[3]);
        tree->Branch("minChisq", &minChisq);
	tree->Branch("dpsi1", &dpsi1);
	tree->Branch("dpsi2", &dpsi2);
	tree->Branch("dpsi3", &dpsi3);

	
	//prepares the final output stream
	ofstream f("../EventOutputs/EventResults.txt");
	
	//iterates minimizes and outputs values over N
	for(int i=0; i<N; i++){
/*for(int i=1; i<=3; i++){
	PrintVector(h->evtP[i].v);
	PrintVector(h->evtP_actual[i].v);
}*/
		evnum=i;
		//reads the next event
		h->populateParticles(h->f1,h->factory1,h->evtP);
		h->populateParticles(h->f2,h->factory2,h->evtP_actual);

		x1g = h->evtP_actual[1].x_m;
                x2g = h->evtP_actual[2].x_m;
                x3g = h->evtP_actual[3].x_m;

		//calculations to be retrived from the minhelper instance from particle structs for use in MINUIT
 		x1m=h->evtP[1].x_m;
		 variance1=mathUtility::getVariance(h->evtP[1]);
		 theta1=h->evtP[1].theta;

		 x2m=h->evtP[2].x_m;
		 variance2=mathUtility::getVariance(h->evtP[2]);
		 theta2=h->evtP[2].theta;

		 x3m=h->evtP[3].x_m;
		 variance3=mathUtility::getVariance(h->evtP[3]);

		 z12=mathUtility::getZ12(h->evtP[1],h->evtP[2]);
		 z23=mathUtility::getZ23(h->evtP[2],h->evtP[3]);
		 z13=mathUtility::getZ13(h->evtP[1],h->evtP[3]);

		//the actual event is read and the minimization is set around the already known min values for x1 x2 so MINUIT doesn't have to work too hard
		 centerX1=h->evtP[1].x_m;
		 centerX2=h->evtP[2].x_m;
		//range for MINUIT to scan for minimum center+/-offset
		offset1=10*sqrt(variance1);
 		offset2=10*sqrt(variance2);
	

		//initializes some MINUIT parameters
		const char * minName = "Minuit2";
		const char *algoName = "DalitzElimination";
		//instance of photon helper to gmathUtility::et the photon fit energy and variance
		PhotonHelper* p = new PhotonHelper(h->evtP);
		
	
		
				//first iteration to initialize a new and old value
				NumericalMinimization(minName,algoName);
			
				x3FitOld = p->getX3Fit(x1Fit,x2Fit);
				std::cout << "x3 Fit value "<< x3FitOld<<std::endl;
				x1FitOld = x1Fit;
				x2FitOld = x2Fit;
			
				variance1 = mathUtility::getVariance(-11,x1FitOld,h->evtP[1].theta);
				variance2 = mathUtility::getVariance(11,x2FitOld,h->evtP[2].theta);
				variance3 = mathUtility::getVariance(22,x3FitOld,-1);	
			
				offset1=10*sqrt(variance1);
 				offset2=10*sqrt(variance2);
				//get a new fit with the variances from the first minimization
				NumericalMinimization(minName,algoName);
				x3Fit = p->getX3Fit(x1Fit,x2Fit);
				std::cout<<"x3 Fit Value "<< x3Fit <<std::endl;
				
			//iterate and compare the fits from each minimization
			while( abs(x3Fit- x3FitOld) > 1e-6 ){
			
				x3FitOld = p->getX3Fit(x1Fit,x2Fit);
				x1FitOld = x1Fit;
				x2FitOld = x2Fit;
				
				variance1 = mathUtility::getVariance(-11,x1FitOld,h->evtP[1].theta);
				variance2 = mathUtility::getVariance(11,x2FitOld,h->evtP[2].theta);
				variance3 = mathUtility::getVariance(22,x3FitOld,-1);
				
				offset1=10*sqrt(variance1);
 				offset2=10*sqrt(variance2);
				NumericalMinimization(minName,algoName);
				x3Fit = p->getX3Fit(x1Fit,x2Fit);
				std::cout<<"x3 Fit Value "<< x3Fit <<std::endl;

			}
			//set x3 fit variance globally, so the calls later are less cluttered
			x3FitVariance = p->getX3Variance(x1Fit,x2Fit,x1FitVariance,x2FitVariance,x12covariance);

		//throw the final data for the pi0 energies onto a histogram
		hESum->Fill(h->evtP[1].v.E()+h->evtP[2].v.E()+h->evtP[3].v.E());
		hESumConstrained->Fill(h->evtP[1].v.E()+h->evtP[2].v.E()+mathUtility::getX3constrained(h->evtP[1].x_m,h->evtP[2].x_m,h->evtP[1],h->evtP[2],h->evtP[3]));

		//calculate new positron and photon energies based on the fit values from the minimization
		TLorentzVector vfit1,vfit2;
		vfit1.SetXYZM((1/x1Fit)*cos(mathUtility::getPhi(h->evtP[1].v)),
			(1/x1Fit)*sin(mathUtility::getPhi(h->evtP[1].v)), 
			(1/x1Fit)*(1/tan(mathUtility::getTheta(h->evtP[1].v))), 
			h->evtP[1].v.M());
	
		vfit2.SetXYZM((1/x2Fit)*cos(mathUtility::getPhi(h->evtP[2].v)), 
			(1/x2Fit)*sin(mathUtility::getPhi(h->evtP[2].v)), 
			(1/x2Fit)*(1/tan(mathUtility::getTheta(h->evtP[2].v))), 
			h->evtP[2].v.M());

		hEfitSum->Fill(vfit1.E()+vfit2.E()+x3Fit);

		EpiPull = ( (vfit1.E()+vfit2.E()+x3Fit) - (h->evtP[1].v.E()+h->evtP[2].v.E()+h->evtP[3].v.E()) ) / 
					sqrt( (variance1+variance2+variance3) - (x1FitVariance+x2FitVariance+x3FitVariance) );
		hEpiPull->Fill(EpiPull);

		//outputs all data to results file
		// data order x_im sigma_x_im x_ireal x_ifit sigma_x_ifit minimum X^2 value
		f<<setprecision(9);
		f<<"3"<<endl;
		f<<h->evtP[1].x_m<<" "<<sqrt(variance1)<<" "<<h->evtP_actual[1].x_m<<" "<<x1Fit<<" "<<sqrt(x1FitVariance)<<" "<<minChisq<<endl;
		f<<h->evtP[2].x_m<<" "<<sqrt(variance2)<<" "<<h->evtP_actual[2].x_m<<" "<<x2Fit<<" "<<sqrt(x2FitVariance)<<" "<<minChisq<<endl;
		f<<h->evtP[3].x_m<<" "<<sqrt(variance3)<<" "<<h->evtP_actual[3].x_m<<" "<<x3Fit<<" "<<sqrt(x3FitVariance)<<" "<<minChisq<<endl;
		
		for(int i=1; i<=3; i++){
			genVects[i]=h->evtP_actual[i].v;
			mesVects[i]=h->evtP[i].v;
		}
		//populate dpsi
		
		dpsi1 = mathUtility::safeAcos( mathUtility::getCosTheta(genVects[1],mesVects[1]) );
		dpsi2 = mathUtility::safeAcos( mathUtility::getCosTheta(genVects[2],mesVects[2]) );
		dpsi3 = mathUtility::safeAcos( mathUtility::getCosTheta(genVects[3],mesVects[3]) );
		
		cout<<setprecision(9);
		if(dpsi3>= 3.14159/2){ 
			cout<<"dpsi: "<< dpsi3<<endl;
			PrintVector(genVects[3]);
			PrintVector(mesVects[3]);
			} 
		hdpsix1->Fill(dpsi1);
		hdpsix2->Fill(dpsi2);
		hdpsix3->Fill(dpsi3);
		
		
		tree->Fill();
	}
	froot->Write();
}
#endif
 
