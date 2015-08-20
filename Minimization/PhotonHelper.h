#ifndef _PHOTON_HELPER
#define _PHOTON_HELPER

#include "../vectorFactory/vectorFactory.h"
#include "../mathUtility/mathUtility.h"
#include <cmath>


//! This Class is a helper class for the MINUIT minimization implementation.  Its function is to calculate the fit energy and uncertainty for the photon from the values generated by the X^2 minimization in MINUIT.
class PhotonHelper{
	public:
	/*! Local copy of the smeared event, should be populated from the minimization implementation and is the same smeared event from the minimization that is read by MinHelper */ 
	vectorFactory::ParticleParameters* evtP;

	/*! Constructor that takes in the smeared event and stores it locally for use in fit calculations */
	PhotonHelper(vectorFactory::ParticleParameters* event);
	//! Function that passes through the minimum x1 and x2 curvatures values for positron and electron to a copy the DalitzChiSq class which calculates the constrained energy for the photon */
	/*!
		\param x1 the curvature of the positron from the X^2 Numerical Minimization
		\param x2 the curvature of the electron from the X^2 Numerical Minimization
		\return The fit energy for the photon based on the X^2 minimization
	*/
	double getX3Fit(double x1, double x2);
	//! Function that uses paramaters of the Covariance matrix in Numerical Minimization to compute the uncertainty on the fitted energy for the photon
	/*!
		\param x1 the curvature of the positron from the X^2 Numerical Minimization
		\param x2 the curvature of the electron from the X^2 Numerical Minimization
		\param x1Variance the variance of the positron from the X^2 Numerical Minimization, retrieved from the covariance matrix
		\param x2Variance the variance of the electron from the X^2 Numerical Minimization, retrieved from the covariance matrix
		\param x12Variance the covariance of the electron and positron curvatures from the X^2 Numerical Minimization, retrieved from the covariance matrix and is less than 0 because the curvatures are anticorrelated
		\return The variance on the photon energy
	*/
	double getX3Variance(double x1, double x2, double x1Variance, double x2Variance, double x12covariance);

	private:
		/*! The mass of the pion */
		double M;
		/*! The mass of the electron/positron */
		double m_e;
		//! Function that assists in calculating the uncertainty on the fit energy.  Computes the partial derivative of the energy fit equation (with mass constraint) with respect to x1
		/*!
			\param x1 the curvature of the positron from the X^2 Numerical Minimization
			\param x2 the curvature of the electron from the X^2 Numerical Minimization
			\return the partial derivative with respect to x1, evaluated with x1 and x2 
		*/
		double getPartialX1(double x1, double x2);
		//! Function that assist in calculating the uncertainty on the fit energy.  Computes the parital derivative of the energy fit equation (with mass constaraint) with respect to x2
				/*!
			\param x1 the curvature of the positron from the X^2 Numerical Minimization
			\param x2 the curvature of the electron from the X^2 Numerical Minimization
			\return the partial derivative with respect to x2, evaluated with x1 and x2 
		*/
		double getPartialX2(double x1, double x2);
		
		
		
	
};
#endif
