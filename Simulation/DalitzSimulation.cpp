#include "../vectorFactory/vectorFactory.h"
#include "../Smearing/DalitzSmear.h"
#include "../Chisq/DalitzChiSq.h"
#include "TLorentzVector.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include "TFile.h"
#include <string>
#include "TH1D.h"
#include "../mathUtility/mathUtility.h"
#include <iomanip>
//debugging single vector printing method
void printVector(TLorentzVector v){
			cout<<v.Px()<<" "<<v.Py()<<" "<<v.Pz()<<" "<<v.E()<<" "<<v.M()<<endl;
	}
//method that updates the calculations in a particle struct after it is smeared
void updateObject(vectorFactory::ParticleParameters& particle){
	vectorFactory* v = new vectorFactory();
	v->makeCalculations(particle);
}
//debugging vector printing method
void printVectors(vectorFactory::ParticleParameters* arrPtr){
		for(int i=1; i<=3; i++){
			cout<<arrPtr[i].pID<<" "<<arrPtr[i].v.Px()<<" "<<arrPtr[i].v.Py()<<" "<<arrPtr[i].v.Pz()<<" "<<arrPtr[i].v.M()<<endl;
	cout<<arrPtr[i].pt<<" "<<arrPtr[i].theta<<" "<<arrPtr[i].x_m<<endl;
		}
	}
int main(int argc, char *argv[]){
	
	double N = atoi(argv[1]);
	double M = atof(argv[2]);
	double m_e = atof(argv[3]);
	double initial_p = atof(argv[4]);
	int chiCont = atoi(argv[5]);
	double scaleParameterNP = atof(argv[6]);
	int seed = atoi(argv[7]);
	


	//In order to histogram the minimalist X^2 values later on we need a TH1D here, it can not be done in the DalitzChiSq class because TH1D can not be accessed as a class member.
	TFile *minimalistFile = new TFile("EventOutputs/DalitzContour.root", "UPDATE");
	TH1D* hminchisq = new TH1D("hminchisq", "Minimalist #chi^{2} ;#chi^{2};Frequency of #chi^{2}",100,0,10);

			//distribution checks for smearing
			TH1D* hepGaus = new TH1D("hepGaus", "x1 Smearing Check; Events per bin;x1m-x1true/sigma_x1m", 100, -5,5);
			TH1D* hemGaus = new TH1D("hemGaus", "x2 Smearing Check; Events per bin;x2m-x2true/sigma_x2m", 100, -5,5); 
			TH1D* hgmGaus = new TH1D("hgmGaus", "x3 Smearing Check; Events per bin;x3m-x3true/sigma_x3m", 100, -5,5);
			TH1D* hE3constrained = new TH1D("hE3constrained", "X3 with pi0 mass constraint; Events per bin; Energy",100,0,10);
			TH1D* hmassconstrained = new TH1D("hmassconstrained", "reconstruction pion mass with x3 constrained; events per bin; mass",100,0.13, 0.14);
	
	
	//print to screen the arguments input to the simulation
	cout<<"ARGS "<<N<<" "<<M<<" "<<m_e<<" "<<initial_p<<" "<<seed<<endl;

	//each time the simulation is run the .hepevt files are wiped for reuse
	ofstream cleaner("EventOutputs/DalitzActualVectors.hepevt");
	cleaner<<"";

	cleaner.close();
	ofstream cleaner2("EventOutputs/DalitzSmearVectors.hepevt");
	cleaner2<<"";
	cleaner2.close();
	//add a cleaner and stream for CM vectors
	ofstream cleaner3("EventOutputs/DalitzCMVectors.hepevt");
	cleaner3<<"";
	cleaner3.close();

	
	ofstream f1("EventOutputs/DalitzActualVectors.hepevt");
	ofstream f2("EventOutputs/DalitzSmearVectors.hepevt");

	//begin simulation loop
	//instansiate vectorfactory to create an event for the simulation and give it the filepath for CM events to pass through to event generator
	vectorFactory* generator =  new vectorFactory(2*m_e/M, "EventOutputs/DalitzCMVectors.hepevt",seed);
	//local particle struct array for the event to be generated
	vectorFactory::ParticleParameters* evtP;
	vectorFactory::ParticleParameters* evtP_actual;
	//single instance of chisq to generate contour if necessary, and histogram the minimalist values
	DalitzChiSq* chisq = new DalitzChiSq();
	//create and smear N events
	for(int i=0; i<N;i++){

		
		//makes an event for this iteration
		
		generator->createEvent(M,m_e,initial_p);
		
		//sets the particle struct array created in event generator to our local evtP pointer
		 evtP = generator->arrPtr;
		 generator->copyObject(evtP,evtP_actual);
		
	
		//write the unsmeared event to the unsmeared event file
		f1<<setprecision(9);
		f1<<3<<endl;
		for(int i=1; i<=3; i++){
			f1<< 1 <<" "<< evtP[i].pID << " " << 0 << " "<< 0 << " " << 0 << " "<< 0 << " ";
			//printVector(evtP[i].v);
			f1<<evtP[i].v.Px()<<" "<<evtP[i].v.Py()<<" "<<evtP[i].v.Pz()<<" "<<evtP[i].v.E()<<" "<<evtP[i].v.M()<<endl;
		}

		//instance of Smearing class to smear the raw event
		DalitzSmear* corruptor = new DalitzSmear();
		//if the directional smearing is not using the default parameters, set them according to input args
		if(scaleParameterNP != 0){ 
			corruptor->setScaleParameterNP(scaleParameterNP);
		}
		//iterate over the particles in the event array, smear each one and then update the particle struct parameters

		//temp lorentz vector to throw out bad direction smears where somehow the dot product between generator and detector 4 vectors <0
		TLorentzVector temp;
		TLorentzVector magTemp;
		for(int i=1;i<=3;i++){

	
			corruptor->setpid(evtP[i].pID);

			
			magTemp = corruptor->SmearVector(evtP[i].v);
			while(mathUtility::getCosTheta(magTemp,evtP[i].v)<0.0){
			//if cos(theta) between the generator and magnitude smeared vector is <0 then the smearing created a anitparallel vector
			//that is, it smeared the energy into the negative range and reversed the direction of the photon
			//so, try smearing it until it gets a good value that doesnt reverse the direction
			
				magTemp = corruptor->SmearVector(evtP[i].v);
				
			}
			evtP[i].v = magTemp;
			
			updateObject(evtP[i]);
			//second smear direction

			//testing infastructure to try to not get negative dot products i.e. psi >= pi/2
			cout<<setprecision(15);
			temp = corruptor->smearDirection(evtP[i].v);
			/*while(mathUtility::getCosTheta(evtP_actual[i].v,temp)<0.0){
				//somehow it didnt work try again, evtP[i] is still generator level at this point
				cout<< "cos theta : "<< mathUtility::getCosTheta(evtP_actual[i].v,temp)<<endl;
				temp = corruptor->smearDirection(evtP[i].v);
				cout<< " new cos theta : "<< mathUtility::getCosTheta(evtP_actual[i].v,temp)<<endl;
			}*/
			// set the good smeared direction vector into the particle object
			evtP[i].v=temp;

	
			updateObject(evtP[i]);
			//extra check to make sure everything went according to plan
			cout<<setprecision(15);
			if(mathUtility::getCosTheta(evtP_actual[i].v,evtP[i].v)<0.0){
				cout<<"BAD DOT PRODUCT    "<<mathUtility::getCosTheta(evtP_actual[3].v,evtP[3].v)<<endl;
				printVector(evtP[i].v);
				printVector(evtP_actual[i].v);
				printVector(magTemp);
			}
			
		}
	
		
		//printVectors(evtP);
		//printVectors(evtP_actual);

		//checking to see if smears are done properly
			 hepGaus->Fill((evtP[1].x_m - evtP_actual[1].x_m)/sqrt(mathUtility::getVariance(evtP[1])));
			 hemGaus->Fill((evtP[2].x_m - evtP_actual[2].x_m)/sqrt(mathUtility::getVariance(evtP[2])));
			 hgmGaus->Fill((evtP[3].x_m - evtP_actual[3].x_m)/sqrt(mathUtility::getVariance(evtP_actual[3])));
			
		//checking to see if energy constraint looks good
		 hE3constrained->Fill(mathUtility::getX3constrained(evtP[1].x_m,evtP[2].x_m,evtP[1],evtP[2],evtP[3]));


		//create histogram of the minimalist X^2
		hminchisq->Fill(chisq->getMinimalistChiSq(evtP[1], evtP[2], evtP[3]));
		//option to create a contour map of the X^2 function recommend only creating the map with small N's e.g. n=1 event
		if(chiCont){
			chisq->generateContour(evtP_actual[1], evtP_actual[2], evtP[1], evtP[2], evtP[3],400, .001);
		}

		//write the smeared event to the smeared event file
		f2<<setprecision(9);
		f2<<3<<endl;
		for(int i=1; i<=3; i++){
			f2<< 1 <<" "<< evtP[i].pID << " " << 0 << " "<< 0 << " " << 0 << " "<< 0 << " ";
			//printVector(evtP[i].v);
			f2<<evtP[i].v.Px()<<" "<<evtP[i].v.Py()<<" "<<evtP[i].v.Pz()<<" "<<evtP[i].v.E()<<" "<<evtP[i].v.M()<<endl;
		}
		

	}
	minimalistFile->Write();

	
	
}
