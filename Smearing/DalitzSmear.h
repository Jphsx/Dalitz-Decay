#ifndef _DALITZ_SMEAR
#define _DALITZ_SMEAR
#include "TLorentzVector.h"
#include "../mathUtility/mathUtility.h"
#include <cmath>
#include "TMath.h"
#include "TRandom.h"
#include "TRandom1.h"
#include <iostream>
#include "TVector3.h"
//#include <math.h>
#include <iomanip>


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

		//static const double scaleParameterCP = 1e-6;
		//static const double scaleParameterNP = 1e-4;
	
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

	//! Method that initiates the sequence of private calls that create a smeared vector v which is the direction smearing of unit vector u contained in input TlorentzVector v
	/*!
		\param v the input four momenta to have direction smeared
		\return direction smeared 4 vector
	*/
	TLorentzVector smearDirection(TLorentzVector v);

	//! method to set the scale parameter for neutral particle i.e. photon, default value is 1e-4
	void setScaleParameterNP(double sigma);
	//! method to set the scale parameter for charged particle i.e. positron/electron, default value is 1e-6
	void setScaleParameterCP(double sigma);
		
		
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



		/////////begin direction smearing///////////////

		//! Gives the unit vector from the associated input 3 vector
		/*!
			\param i i component of (i,j,k)
			\param j j component of (i,j,k)
			\param k k component of (i,j,k)
			\return the unit vector (i,j,k)
		*/
		TVector3 getUnitVector(double i, double j, double k);
		//! Gives the magnitude vector of the input 3 vector
		/*!
			\param i i component of (i,j,k)
			\param j j component of (i,j,k)
			\param k k component of (i,j,k)
			\return vector magnitude
		*/
		double getVectorMagnitude(double i, double j, double k);
		//! calculates the omega_t unit vector, which is orthogonal to the original input vector to be smeared. The omega_t vector dictates the direction on the plane normal to the orignal vector will be smeared.
		/*!
			\param v_1 the quantity returned from getV_1 
			\param v_n a vector that is orthogonal to v_1 and the original unit vector u
			\return the smearing direction unit 3 vector omega_t
		*/
		TVector3 getOmega_T(TVector3 v_1, TVector3 v_n);
		//! finds a vector v1 orthogonal to the original input vector u, and is intended to be a component of omega_t
		/*!
			\param uVector the original, to be smeared, unit vector
			\return the omega_t component v_1 3 vector
		*/
		TVector3 getV_1(TVector3 uVector);
		//! calculates the dot product between to vectors v1 . v2
		double getScalarProduct(TVector3 v1, TVector3 v2);
		//! calculates the cross product between to vectors v1 x v2 
		TVector3 getVectorProduct(TVector3 v1, TVector3 v2);
		//! calculates the angular deviaton dpsi for vector v based on a random deviate drawn from a Rayleigh distribution that depends on scale parameter sigma. The scale parameters selected depend on the particle's PID owned by the class at the time of smearing. i.e. the scale parameters for charged particles and neutral particles should be different.
		/*!
			\param v The four momenta which direction is to be smeared.
			\return angular deviation dpsi
		*/
		double getDpsi(TLorentzVector v);

		//! The scale parameter sigma for a charged particle, for drawing random deviates from rayleigh distribution 
		double scaleParameterCP;
		//! The scale parameter sigma for a neutral particle, for drawing random deviates from rayleigh distribution
		double scaleParameterNP;

		
		
};
#endif
