#ifndef _VECTOR_FACTORY
#define _VECTOR_FACTORY
#include <iostream>
#include "TLorentzVector.h"
#include "TFile.h"
#include <cmath>
#include <string>
#include <fstream>
#include "TMath.h"


using namespace std;
//reads 4vector from file, eventually should make them as well
class vectorFactory{
	public:
	void readVector(const char* path); 


	struct ParticleParameters{
		TLorentzVector v;
		double theta;
		double pt;
		double x_m;
		int pID;	
		
		};
	//indices should match particle numbers
	ParticleParameters* arrPtr;

	ParticleParameters* getParticles();

	//struct population functions
	void setTLVector(int pID,double px,double py,double pz, double m);
	void makeCalculations(ParticleParameters& particle);
	

	//debugging
	void printVectors(){
		for(int i=1; i<=3; i++){
			cout<<arrPtr[i].pID<<" "<<arrPtr[i].v.Px()<<" "<<arrPtr[i].v.Py()<<" "<<arrPtr[i].v.Pz()<<" "<<arrPtr[i].v.M()<<endl;

			cout<<arrPtr[i].pt<<" "<<arrPtr[i].theta<<" "<<arrPtr[i].x_m<<endl;
		}
		
	}
	

};
#endif
