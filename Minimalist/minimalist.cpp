#ifndef MINIMALIST_
#define MINIMALIST_

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


using namespace std;

int main(int argc, char* argv[]){

int N  = atoi(argv[1]);
double M = atof(argv[2]);
double m_e = atof(argv[3]);
double initP = atof(argv[4]);

TFile *froot = new TFile("../EventOutputs/ResultHistos_Minimalist.root", "UPDATE");

TH1D* hESum = new TH1D("hESum", "Measured pi0 Energy E1+E2+E3; GEV; Event Per Bin",100,initP-0.5*initP,initP+0.5*initP);
TH1D* hEminSum = new TH1D("hEminSum", "Minimalist pi0 Energy E1+E2+ E3(x1,x2); GEV; Event Per Bin", 500, initP-3.5*initP,initP+3.5*initP);
TH1D* hEpiPull = new TH1D("hEpiPull", "pull distribution of pi0 energies; pull value; event per bin",100,-0.5*initP,0.5*initP);
TH1D* hmass_m = new TH1D("hmass_m", "pi0 measured mass resolution; Mass; event per bin", 100, M-M,M+M);
TH1D* hmass_min = new TH1D("hmass_min", "pi0 minimalist mass resolution; Mass; event per bin", 100, M-M, M+M);
//DEBUGGING
ofstream debug("../EventOutputs/Event_min_debug.txt");




MinHelper* h = new MinHelper("../EventOutputs/DalitzSmearVectors.hepevt","../EventOutputs/DalitzActualVectors.hepevt");




ofstream cleaner("../EventOutputs/EventResults_Minimalist.txt");
cleaner<<"";
cleaner.close();

ofstream f("../EventOutputs/EventResults_Minimalist.txt");

	double chisq;
	double E3constrained;
	double ESum_m;
	TLorentzVector massV,vtest;
	double mass;

	for(int i=0; i<N; i++){

		h->populateParticles(h->f1,h->factory1,h->evtP);
		h->populateParticles(h->f2,h->factory2,h->evtP_actual);

				debug<<setprecision(9);
				debug<<"event# "<<i<<endl;
			for(int j=1; j<=3; j++){
				debug<<"pid"<<h->evtP[j].pID <<" "
				<<"E: "<<h->evtP[j].v.E() << " "
				<<"M: "<< h->evtP[j].v.M()<< " "
				<<"P: "<<h->evtP[j].v.P()<< " "
				<<"x_m: "<< h->evtP[j].x_m << " "
				<<"theta: "<<h->evtP[j].theta<<endl<<endl;
			}
			 E3constrained = mathUtility::getX3constrained(h->evtP[1].x_m,h->evtP[2].x_m,h->evtP[1],h->evtP[2],h->evtP[3]);

			//debug constraint equations
			 debug <<"E3C: "<< E3constrained <<" "
			<<"E1+E2+E3C: "<< h->evtP[1].v.E() + h->evtP[2].v.E() + E3constrained 
			<< "E1+E2+E3: " << h->evtP[1].v.E() + h->evtP[2].v.E() + h->evtP[3].v.E()<< endl<<endl;;
			//
		//currently uncertainty is just from E3m detector level, the uncertainty may need to be propogated from x1, x2
		chisq = pow(E3constrained- h->evtP[3].x_m,2)/mathUtility::getVariance(h->evtP[3]);

//write the new results file
// data order x_im sigma_x_im x_ireal x_imin sigma_x_imin minimum X^2 value
//NOTE: variance on x1x2 are the same so may return infinity in pull dist, same with x3, need to propagate error on x3

		//here is the error propagation on x3 from x1 x2
		double varianceX3;
		PhotonHelper* dx3 = new PhotonHelper(h->evtP);
		varianceX3=dx3->getX3Variance(h->evtP[1].x_m, h->evtP[1].x_m , mathUtility::getVariance(h->evtP[1]), mathUtility::getVariance(h->evtP[2]) , 0);
		

		f<<setprecision(9);
		f<<"3"<<endl;
		f<<h->evtP[1].x_m<<" "<<sqrt(mathUtility::getVariance(h->evtP[1]))<<" "<<h->evtP_actual[1].x_m<<" "<<h->evtP[1].x_m<<" "<<sqrt(mathUtility::getVariance(h->evtP[1]))<<" "<<chisq<<endl;
		f<<h->evtP[2].x_m<<" "<<sqrt(mathUtility::getVariance(h->evtP[2]))<<" "<<h->evtP_actual[2].x_m<<" "<<h->evtP[1].x_m<<" "<<sqrt(mathUtility::getVariance(h->evtP[2]))<<" "<<chisq<<endl;
		f<<h->evtP[3].x_m<<" "<<sqrt(mathUtility::getVariance(h->evtP[3]))<<" "<<h->evtP_actual[3].x_m<<" "<<E3constrained<<" "<<sqrt(varianceX3)<<" "<<chisq<<endl;

	
		ESum_m = h->evtP[1].v.E() + h->evtP[2].v.E() + h->evtP[3].v.E();
		hESum->Fill(ESum_m);
		hEminSum->Fill(h->evtP[1].v.E() + h->evtP[2].v.E() + E3constrained);
		
		massV = h->evtP[1].v+h->evtP[2].v+h->evtP[3].v;
		hmass_m->Fill(massV.M());
	
		//vtest.SetPxPyPzE(E3constrained,0,0,E3constrained);
		//massV = h->evtP[1].v +h->evtP[2].v + vtest;
		mass =  sqrt( (pow(h->evtP[1].v.E()+h->evtP[2].v.E()+E3constrained,2)) - (pow(h->evtP[1].v.P()+h->evtP[2].v.P()+E3constrained,2)) );
		//hmass_min->Fill(mass);
		hmass_min->Fill(mass);	
	}

	//hEpiPull->Fill(  (ESum_m - (h->evtP[1].x_m + h->evtP[2].x_m +E3constrained))/(sqrt( // what is the epi0 variance??
	
	

	froot->Write();


}
#endif
