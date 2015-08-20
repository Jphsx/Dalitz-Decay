#ifndef _DALITZ_MATH_
#define _DALITZ_MATH_
#include "TLorentzVector.h"
#include "../vectorFactory/vectorFactory.h"
#include <cmath>

typedef vectorFactory::ParticleParameters pObject;
//! Namespace that combines all the common Dalitz Math Utilities overloaded to support multiple class usage e.g. finding variances or angular quantities of vectors
namespace mathUtility {


			const double M=0.13497;
			const double m_e=0.000511;
			//! Given a particle struct with a defined particle ID, this function returns the variance of that particle's measured curvature if it is a photon or electron, or the variance of the measured energy of the photon. Based on the assumption that the collider is cylindrical.
			/*!
			\param p a valid particle struct with a defined particle ID
			*/
			 double getVariance(pObject p);

			//! Returns the variance of the input vector with a previously class defined particle ID.  Same variance definition as the DalitzChiSq class. Used solely for finding the sigmas for gaussian number generation.
			/*!
			\param v The input "v_reg" 4 vector
			*/
			 double getVariance(TLorentzVector v, int pid);
			
			//! Returns the variance of particle based on the arguments given, 11,-11 will give a variance electron positron curvature and 22 will get the variance on the photon energy.  The variances are calculated based on the x value given. This function is overloaded to get uncertainties on measured values based on x values that come from fitting
			/*!
			\param pID the particle ID
			\param xm the measured x value (xfit is the intended use)
			\theta the polar angle of the particle with respect to the z axis in the accelerator
			
			*/
			 double getVariance(int pID, double xm, double theta);

			 //! Evaluates a function based on the opening angle between the positron and electron, used for evaluating the Energy of the photon in the event in terms of x1, x2 and the pion mass constraint
			/*!
			\param p1 the particle struct for the positron
			\param p2 the particle struct for the electron
			*/
			 double getZ12(pObject p1, pObject p2);
			 
			 //! Evaluates a function based on the opening angle between the photon and electron, used for evaluating the Energy of the photon in the event in terms of x1, x2 and the pion mass constraint
			/*!
			\param p2 the particle struct for the electron
			\param p3 the particle struct for the photon
			*/
			 double getZ23(pObject p2, pObject p3);
			
			 //! Evaluates a function based on the opening angle between the positron and photon, used for evaluating the Energy of the photon in the event in terms of x1, x2 and the pion mass constraint
			/*!
			\param p1 the particle struct for the positron
			\param p3 the particle struct for the photon
			*/
			 double getZ13(pObject p1, pObject p3);

			 //! A function that calculates the opening angle between two vectors v1 and v2
			/*!
			\param v1 a TLorentzVector 4vector
			\param v2 a TLorentzVector 4vector
			\return opening angle in radians between the two vectors
			*/
			 double getCosTheta(TLorentzVector v1, TLorentzVector v2);


			
			//! Function that gets the polar angle of the pz/p
			/*!
			\param v The input 4 vector
			*/
			 double getTheta(TLorentzVector v);

			//! Calculates the scalar magnitude of the input vector v
			/*!
			\param v The input 4 vector
			*/
			 double get_p(TLorentzVector v);

			//! Calculates the azimuthal angle of the input vector v
			/*!
			\param v The input 4 vector
			*/
			 double getPhi(TLorentzVector v);

			//! Function that gets the scalar magnitude of the x and y components of an input 4vector
			/*!
			\param v The input 4 vector
			*/
			 double getPt(TLorentzVector v);

			//! A function that evaluates the energy of the photon in the event in terms of the positron and electron curvatures, angles between the particles and uses the pion mass constraint.
			/*!
			\param x1 the curvature of the positron
			\param x2 the curvature of the electron
			\param p1 the particle struct of the positron
			\param p2 the particle struct of the electron
			\param p3 the particle struct of the photon
			\return the photon constrained energy
			*/
			 double getX3constrained(double x1,double x2, pObject p1, pObject p2, pObject p3);
			
};
#endif
