#include <iostream>
#include "TRandom1.h"
#include <cmath>
#ifndef MC_GENERATOR_
#define MC_GENERATOR_
using namespace std;
//! Generates x,y values that govern the properties of the dielectric mass, such that x is the invariant mass of the dielectric system and y is the Energy partition in the CM frame
class XYgenerator{
	
	public:
		//! generates a valid X,Y pair based on the globally defined boundary conditions set by the constructor, and uses the test probability function to validate the values generated.
		/*!
			\return a pointer to a double array containing x and index 0 and y at index 1
		*/
		double* generateXYpair();

		//! Constructor that initializes the random number generator and sets the global bottom boundary for the x value
		/*!
			\param bound or nu=2m_e/M and defines the boundary such that nu^2<x<1
		*/
		XYgenerator(double bound);

		//! Constructor that utilizes an external random number generator and sets the global bottom boundary for the x value
		/*!
			\param bound or nu=2m_e/M and defines the boundary such that nu^2<x<1
		*/
		XYgenerator(double bound, TRandom1* rng);
		
	
		//! Calculates and provides the beta values intended for y such that -beta<y<beta based on a pregenerated x value
		/*!
			\param x the invariant dielectric mass
			\return the positive upper bound for y
		*/
		double getBeta(double x);
		//! Calculates the probability of the differential decay rate given the parameters x and y
		/*!
			\param x the invariant dielectric mass
			\param y the energy partition in the CM frame
			\return the value of P(x,y)
		*/
		double getProbability(double x, double y);
		//! The Monte Carlo rejection algorthim returns true for a valid x,y set whose probability is greater than a random number generated between 0 and WMAX
		/*!
			\param p the probability evaluated from the get probability function
			\param x the invariant dielectric mass
			\param y the energy partition in the CM frame
			\return true or false based on whether the pair of points is valid or rejected, respectively
		*/
		bool testProbability(double p, double x, double y);
		/*! The approximate maximum of the diferential decay rate (probability density function) */
		double WMAX; 
		/*! nu=2m_e/M and defines the boundary of the dielectric invariant mass such that nu^2<x<1 */
		double nu; 
		/*! Local instance of random number generator that is used for creating x,y values and test points for the rejection algorthim */
		TRandom1 *RNG; 
		
};
#endif
