//#include "vectorFactory.h"
#ifndef _EVENT_GEN
#define _EVENT_GEN
//#define _USE_MATH_DEFINES
#include "TLorentzVector.h"
#include "TRandom1.h"
#include "../MC_Rejection/XYgenerator.h"
#include <cmath>
#include <math.h>
#include "TVector3.h"
#include <iostream>
#include <fstream>



using namespace std;
//! Creates unique three particle events in the CM frame and Lab frame
class eventGenerator{
	public:
		//! Creates a three particle event and returns the set of particles in a 4 vector array.
		/*!
			\param M the mass of the pion.
			\param m_e the mass of the electron/positron.
			\param initial_p The initial momentum of the system in the lab frame.
			\return A pointer to TLorentzVector array of length N particles+1.
		*/
		TLorentzVector* generateEvent( double M, double m_e, double initial_p);
		//! Constructor to initialize all local class members such as the random number generator and Monte Carlo x,y generator.
		/*!
			\param nu  2m_e/M where nu^2 is the bottom boundary for the parameter x (nu used in the constructor argument of XY generator, this constructor is simply a pass through).
		*/
		eventGenerator(double nu);
		
		//! Constructor to initialize all local class members such as the random number generator and Monte Carlo x,y generator, this constructor is used specifically to supply the class with a filepath output CM events into
		/*!
			\param nu 2m_e/M where nu^2 is the bottom boundary for the parameter x (nu used in the constructor argument of XY generator, this constructor is also a pass through)
			\param CMfilename the string path used for the filestream to output CM events in the generateEvent() method
		*/
		eventGenerator(double nu, const char* CMfilename);

		//! Constructor to initialize all local class members such as the random number generator and Monte Carlo x,y generator, this constructor is used specifically to supply the class with a filepath output CM events into as well as a seed for xy generation and event generation
		/*!
			\param nu 2m_e/M where nu^2 is the bottom boundary for the parameter x (nu used in the constructor argument of XY generator, this constructor is also a pass through)
			\param CMfilename the string path used for the filestream to output CM events in the generateEvent() method
			\param seed the seed to be used for TRandom1 which is the RNG object for xygenerator and event generator
		*/
		eventGenerator(double nu, const char* CMfilename, int seed);
		

	private:
		//! Takes in a 4vector array and rotates the particle system about the x,y and z axis.
		/*!
			\param evtP address of the pointer of a TLorentzVector array.
		*/
		void rotateSystem(TLorentzVector*& evtP);
		//! Boosts the particles in the referenced 4vector array along a boost vector whose magnitude is determined by the initial conditions but has a randomly generated direction.
		/*!
			\param evtP The address of the pointer to a 4 vector array.
			\param initial_p The initial momentum of the system in the lab frame.
			\param M the mass of the pion.
		*/
		void boostToLab(TLorentzVector*& evtP, double initial_p, double M);
		//! Creates a 4vector for a photon in the event, given the conditions provided by the arguments.
		/*!
			\param x the electron/positron invariant mass, created by XYgenerator.
			\param M the mass of the pion.
			\return a TLorentzVector for the photon.
		*/	
		TLorentzVector makePhoton(double x, double M);
		//! Creates a 4vector for an electron in the event, given the conditions provided by the arguments.
		/*!
			\param x the electron/positron invariant mass, created by XY generator.
			\param y Energy partition in the CM frame, created by XY generator.
			\param M the mass of the pion.
			\param m_e the mass of the electron/positron.
			\return a TLorentzVector for the electron.
		*/
		TLorentzVector makeElectron(double x, double y,double M, double m_e);
		//! Creates a 4vector for a positron in the event, given the conditions provided by the arguments.
		/*!
			\param x the electron/positron invariant mass, created by XY generator.
			\param y Energy partition in the CM frame, created by XY generator.
			\param M the mass of the pion.
			\param m_e the mass of the electron/positron.
			\return a TLorentzVector for the positron.
		*/
		TLorentzVector makePositron(double x, double y,double M, double m_e);
		/*! Local copy of the random number generator used for generating the random directions of the particle rotations and random boost direction.*/
		TRandom1* RNG;
		/*! Local copy of XYgenerator to supply particle creation methods with x and y values. */
		XYgenerator* xygen;

		//! Debugging method that prints the contents of a 4vector array.
		void vectorPrinter(TLorentzVector* v);

		/*! Filestream shunted into event generator to output the CM events to .hepevt files */
		ofstream CMfs;
	
		/*! Boolean flag to indicate to the class whether to write CM events to file */
		bool CMwrite;

		//! Helper method to write an 4vector Array in the specified file steam, the PIDs for the standard .hepevt output format are hardcoded into a array that exists only in the scope of this function
		/*!
			\param f1 the filestream to be written to.
			\param evtP the CM rotated 4vector array
		*/
		void writeEvent(ostream& f1, TLorentzVector*& evtP);
		
		
};
#endif
