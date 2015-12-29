
#ifndef _VECTOR_FACTORY
#define _VECTOR_FACTORY
#include <iostream>
#include "TLorentzVector.h"
#include "TFile.h"
#include <cmath>
#include <string>
#include <fstream>
#include "TMath.h"
#include "../eventGenerator/eventGenerator.h"


using namespace std;

//! Reads events from a specified file or calls event generator to create a 3 particle event and packages a raw event into a particle struct array. 
class vectorFactory{
	public:
		
		//! Event particle container
		/*! Contains the particle 4 vector, theta, transverse momentum, measured curvature, and the particle ID. For photons the x_m value is the measured energy and theta/pt are unused.*/
		struct ParticleParameters{
		/*! TLorentz particle 4 vector for a raw unsmeared event */
		TLorentzVector v; 
		/*! Particle's polar angle with respect to z-axis (electron beam direction) *Note* Photon theta is unused in all calculations, so it remains initialized to 0.0.  */
		double theta; 
		/*! The particle's transverse momentum. */
		double pt; 
		/*! The particle's measured curvature (1/pt). */
		double x_m; 
		/*! The particle ID number. */	
		int pID; 
		
		};
		

	/*! Local copy of particle struct array, indices match the standard particle numbers i.e. 1-positron, 2-electron, 3-photon. */
	ParticleParameters* arrPtr;
	/*! Local instance of eventGenerator to get 4 vectors */
	eventGenerator* gen;
	
	//ParticleParameters* getParticles();

	//! (deprecated with N events) A function that reads a specified .hepevt file and calls the struct population functions.
	/*!
		\param path The path the the specified .hepevt file.
	*/
	void readVector(const char* path);
	//! A function that calls the local eventGenerator to create a unique 3 particle event and calls the struct population functions
	/*!
		\param M The mass of the pion.
		\param m_e The mass of the electron/positron.
		\param initial_p The initial momentum of the system in the Lab frame.
	*/ 
	void createEvent(double M, double m_e, double initial_p);
	//! Empty constructor for quick use of struct population functions that do not require event generation
	vectorFactory();
	//! Constructor used for generating N events
	/*!
		\param nu 2m_e/M which the square is used as the bottom boundary for generating x,y values with the Monte Carlo rejection method
	*/
	vectorFactory(double nu);
	//! Constructor used for generating N events and creating a output .hepevt file which contains the CM events
	/*!
		\param nu 2m_e/M which the square is used as the bottom boundary for generating x,y values with the Monte Carlo rejection method
		\param CMfilepath the filepath designated for eventGenerator to output CM events into
	*/
	vectorFactory(double nu, const char* CMfilepath);

	//! Constructor used for generating N events and creating a output .hepevt file which contains the CM events, uses a externally defined seed for eventGenerator
	/*!
		\param nu 2m_e/M which the square is used as the bottom boundary for generating x,y values with the Monte Carlo rejection method
		\param CMfilepath the filepath designated for eventGenerator to output CM events into
		\param seed the seed for the RNG object in eventGenerator and XYgenerator
	*/
	vectorFactory(double nu, const char* CMfilepath, int seed);
	

	//struct population functions
	//! Given the vector information and a particle ID, this helper method will populate the 4 vector for a struct whose index is associated with the particle ID
	/*!
		\param pID The particle ID which determines the index of where the 4 vector will be set.
		\param px The x-component of momentum.
		\param py The y-component of momentum.
		\param pz The z-component of momentum.
		\param m The mass of the particle.
	*/
	void setTLVector(int pID,double px,double py,double pz, double m);
	//! Takes a reference address to a particle struct with a predefined particle ID and sets the struct members at the address (pt, theta, x_m).
	/*!
		\param particle a particle struct passed by reference.
	*/
	void makeCalculations(ParticleParameters& particle);
	//! Given a particle struct array, sets all the particle IDs at the appropriate indices.
	/*!
		\param evtP event particles, a pointer to a particle struct array.
	*/
	void setPID(ParticleParameters* evtP);
	
	//! Method that makes a deep copy of one struct array onto another struct array
	/*!
		\param evtP_a the particle struct array to be copied
		\param evtP_b the particle struct array to become a copy of the first argument 
	*/
	void copyObject(ParticleParameters* evtP_a, ParticleParameters* evtP_b);

	//debugging
	//! Debugging method that prints the contents of the local particle struct array.
	void printVectors(){
		for(int i=1; i<=3; i++){
			cout<<arrPtr[i].pID<<" "<<arrPtr[i].v.Px()<<" "<<arrPtr[i].v.Py()<<" "<<arrPtr[i].v.Pz()<<" "<<arrPtr[i].v.M()<<endl;

			cout<<arrPtr[i].pt<<" "<<arrPtr[i].theta<<" "<<arrPtr[i].x_m<<endl;
		}
		
	}

		
	

};
#endif
