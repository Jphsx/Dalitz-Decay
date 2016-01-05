#ifndef _MINIMALIST2_
#define _MINIMALIST2_

#include "../Minimization/MinHelper.h"
#include "../Minimization/PhotonHelper.h"
#include <iostream>
#include <cmath>
#include "TLorentzVector.h"
#include <fstream>
#include "TFile.h"
#include "TH1D.h"
#include "../mathUtility/mathUtility.h"
#include <iomanip>
#include "min2utility.h"

using namespace std;



int main(int argc, char* argv[]){
int N  = atoi(argv[1]);
double M = atof(argv[2]);
double m_e = atof(argv[3]);
double initP = atof(argv[4]);
//double scalePcp = atof(argv[5]);

TFile *froot = new TFile("../EventOutputs/ResultHistos_Minimalist2.root", "UPDATE");

TH1D* hESum = new TH1D("hESum", "Measured pi0 Energy E1+E2+E3; GEV; Event Per Bin",100,initP-0.5*initP,initP+0.5*initP);
TH1D* hEmin2Sum = new TH1D("hEmin2Sum", "Minimalist pi0 Energy E1+E2+ E3(x12); GEV; Event Per Bin", 500, initP-3.5*initP,initP+3.5*initP);
TH1D* hEpiPull = new TH1D("hEpiPull", "pull distribution of pi0 energies; pull value; event per bin",100,-0.5*initP,0.5*initP);

TH1D* hmass_min2 = new TH1D("hmass_min2", "pi0 minimalist2 mass resolution; Mass; event per bin", 100, M-M, M+M);
TH1D* hmass_m = new TH1D("hmass_m", "pi0 measured mass resolution; Mass; event per bin", 100, M-M,M+M);

//DEBUGGING
ofstream debug("../EventOutputs/Event_min2_debug.txt");

MinHelper* h = new MinHelper("../EventOutputs/DalitzSmearVectors.hepevt","../EventOutputs/DalitzActualVectors.hepevt");




ofstream cleaner("../EventOutputs/EventResults_Minimalist2.txt");
cleaner<<"";
cleaner.close();

ofstream f("../EventOutputs/EventResults_Minimalist2.txt");

	double chisq;
	double E3constrainedM2;
	double ESum_m;
	double psi123est;
	double massC,massM;
	//TLorentzVector v12;
	min2utility* min2helper = new min2utility( M );
	for(int i=0; i<N; i++){


		h->populateParticles(h->f1,h->factory1,h->evtP);
		h->populateParticles(h->f2,h->factory2,h->evtP_actual);
	min2helper->setVectors(h->evtP[1].v+h->evtP[2].v, h->evtP[3].v);
		//use bisection on the X^2 derivative, the zero is the optimal opening angle such that X^2 is minimized
		//first set the sigma of opening angle since it is now 1mrad/rootE
		
		min2helper->setSigma_psi12_3(.001/sqrt(h->evtP[3].v.E()));
		psi123est=min2helper->MinimizeMin2();

		E3constrainedM2 = mathUtility::getX3constrainedMin2(h->evtP[1].v+h->evtP[2].v,h->evtP[3].v, psi123est);
		
		chisq = min2helper->getChiSqMin2(M, h->evtP[1].v+h->evtP[2].v, h->evtP[3].v ,  psi123est, mathUtility::safeAcos(mathUtility::getCosTheta(h->evtP[1].v+h->evtP[2].v,h->evtP[3].v)));


		//debugging framework
						debug<<"event# "<<i<<endl;
			for(int j=1; j<=3; j++){
				debug<<"pid"<<h->evtP[j].pID <<" "
				<<"E: "<<h->evtP[j].v.E() << " "
				<<"M: "<< h->evtP[j].v.M()<< " "
				<<"x_m: "<< h->evtP[j].x_m << " "
				<<"theta: "<<h->evtP[j].theta<<endl<<endl;
				
			}
			 

			//debug constraint equations
			debug<<"PSI12,3 Real: "<< mathUtility::safeAcos(mathUtility::getCosTheta(h->evtP_actual[1].v+h->evtP_actual[2].v,h->evtP_actual[3].v))<<" "
			<<"PSI12,3 Meas: "<< mathUtility::safeAcos(mathUtility::getCosTheta(h->evtP[1].v+h->evtP[2].v,h->evtP[3].v))<<" "
			<<"PSI12,3 Adjusted: "<< psi123est<< " CHISQ = "<<chisq<<endl<<endl;
			 debug <<"E3C: "<< E3constrainedM2 <<" "
			<<"E1+E2+E3C: "<< h->evtP[1].v.E() + h->evtP[2].v.E() + E3constrainedM2 
			<< " E1+E2+E3: " << h->evtP[1].v.E() + h->evtP[2].v.E() + h->evtP[3].v.E()<< endl<<endl;;
	//////////////////

	//need to add to photon helper total differential
	//for min2

	double varianceX3=1.0;
		/*PhotonHelper* dx3 = new PhotonHelper(h->evtP);
		varianceX3=dx3->getX3Variance(h->evtP[1].x_m, h->evtP[1].x_m , mathUtility::getVariance(h->evtP[1]), mathUtility::getVariance(h->evtP[2]) , 0); */
		
//update all these outputs because they are from min1

		f<<setprecision(9);
		f<<"3"<<endl;
		f<<h->evtP[1].x_m<<" "<<sqrt(mathUtility::getVariance(h->evtP[1]))<<" "<<h->evtP_actual[1].x_m<<" "<<h->evtP[1].x_m<<" "<<sqrt(mathUtility::getVariance(h->evtP[1]))<<" "<<chisq<<endl;
		f<<h->evtP[2].x_m<<" "<<sqrt(mathUtility::getVariance(h->evtP[2]))<<" "<<h->evtP_actual[2].x_m<<" "<<h->evtP[1].x_m<<" "<<sqrt(mathUtility::getVariance(h->evtP[2]))<<" "<<chisq<<endl;
		f<<h->evtP[3].x_m<<" "<<sqrt(mathUtility::getVariance(h->evtP[3]))<<" "<<h->evtP_actual[3].x_m<<" "<<E3constrainedM2<<" "<<sqrt(varianceX3)<<" "<<chisq<<endl;

	
		ESum_m = h->evtP[1].v.E() + h->evtP[2].v.E() + h->evtP[3].v.E();
		hESum->Fill(ESum_m);
		hEmin2Sum->Fill(h->evtP[1].v.E() + h->evtP[2].v.E() + E3constrainedM2);
	
		massM = mathUtility::getConstrainedMass(h->evtP[1].v , h->evtP[2].v , h->evtP[3].v);		
		hmass_m->Fill(massM);		

		massC = mathUtility::getConstrainedMass(h->evtP[1].v + h->evtP[2].v , E3constrainedM2, psi123est);
		hmass_min2->Fill(massC);

	}
	froot->Write();
		


}
#endif
