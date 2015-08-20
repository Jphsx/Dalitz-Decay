#ifndef _DALITZ_SMEAR
#define _DALITZ_SMEAR
#include "TLorentzVector.h"
#include <cmath>
#include "TMath.h"
#include "TRandom.h"
#include <iostream>


using namespace std;

class DalitzSmear{
	public:
		
		TRandom *RNG;
		TLorentzVector v_reg;
		TLorentzVector v_sme;
		double q;
		int pid;
	
		void setpid(int id);

	//constructor
	DalitzSmear();
	DalitzSmear(TLorentzVector v, int classifier);
	
	//multiple ways to smear a vector

	//if default empty constructor pass in vector&
	//if populated from constructor, just returns local v_sme
	TLorentzVector SmearVector();
	TLorentzVector SmearVector(TLorentzVector v);	

		
	private://not for access by other classes local calculations only
		double getPt(TLorentzVector v);
		double getTheta(TLorentzVector v);
		double get_p(TLorentzVector v); //p^2 = pt^2 + pz^2
		double getPhi(TLorentzVector v);
		double getPtsm(TLorentzVector v);

		double getGauss(TLorentzVector v);
		double getVariance(TLorentzVector v);
		

};
#endif
