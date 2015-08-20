#ifndef _DALITZ_SMEAR
#define _DALITZ_SMEAR
#include "TLorentzVector.h"
#include "../mathUtility/mathUtility.h"
#include <cmath>
#include "TMath.h"
#include "TRandom.h"
#include "TRandom1.h"
#include <iostream>


using namespace std;
//! Class that takes in a 4vector and smears it, each vector both smeared and unsmeared is stored locally, and the class requires the particle id to be seperately supplied to the class (where it is also stored locally).
class DalitzSmear{
	public:
		//! Local copy of random number generator used for generating gaussian values
		TRandom1 *RNG;
		//! Locally stored 4 vector which is unsmeared
		TLorentzVector v_reg;	
		//! Locally stored 4 vector which is the smeared version of v_reg
		TLorentzVector v_sme;
		//! Charge, used only electron(-1) and positron(+1)
		double q;
		//! The local particle ID associated with the current copy of v_reg
		int pid;
	
		//! Method to set the particle ID
		/*!
			\param id the particle ID associated with the 4 vector that is intended to be smeared
		*/
		void setpid(int id);

	//constructor
	//! Empty constructor
	DalitzSmear();
	//! Constructor that sets the 4 vector to be smeared, and the particle ID associated with the 4 vector, highly recommended to use this constructor over the empty version.
	DalitzSmear(TLorentzVector v, int classifier);
	
	//multiple ways to smear a vector

	//if default empty constructor pass in vector&
	//if populated from constructor, just returns local v_sme
	//! Method that smears a 4 vector *NOTE* v_reg and pid global variables must already be set.
	TLorentzVector SmearVector();
	//! Method that smears the argument 4 vector *NOTE* the particle ID must be set before calling this method.
	/*!
		\param v The 4vector smearing target
	*/
	TLorentzVector SmearVector(TLorentzVector v);	

		
	private://not for access by other classes local calculations only
		//! Function that returns the smeared values for Pt, which is used for smearing each component of the v_reg vector
		/*!
			\param v The input "v_reg" 4 vector
		*/
		double getPtsm(TLorentzVector v);
		//! Generates a random gaussian value for the input vector. If the associated particle ID is an electron or positron the gaussian mean is the curvature of the particle with its associated sigma.  If the associated partice ID is a photon the gaussian mean is the measured photon energy with its associated sigma
		/*!
			\param v The input "v_reg" 4 vector
		*/
		double getGauss(TLorentzVector v);
		
		

};
#endif
