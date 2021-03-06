#ifndef _DALITZ_CHI_SQ
#define _DALITZ_CHI_SQ
#include "TLorentzVector.h" 
#include "../vectorFactory/vectorFactory.h"
#include <cmath>
#include "TFile.h"
#include "TH2.h"
#include "../mathUtility/mathUtility.h"


typedef vectorFactory::ParticleParameters pObject;
//! Class that generates elimination X^2 contours and minimalist X^2 values, also used as a calculation helper for other classes to retrieve the variances or opening angles of and between particles. Uses the particle structs from the vectorFactory class as the standard for accessing information about the event and particles
class DalitzChiSq
{
	friend class MinHelper;
	
	public:
		/*! the pion mass */
		double M;
		/*! the electron positron mass */
		double m_e;
		/*! global copy the photon energy computed with the pion mass constraint (unused) */
		double x3constrained;		

		//! Evaluates the X^2(x1,x2,x3) function such that x1 x2 X^2 terms are exactly 0 and the mean x3 is given by the mass constrained energy
		/*!
			\param p1 the particle 1 (positron) struct
			\param p2 the particle 2 (electron) struct
			\param p3 the particle 3 (photon) struct
			\return the X^2(x1,x2,x3) value evaluated such that x1 and x2 have 0 contribution and x3(x1,x2)
		*/
		double getMinimalistChiSq(pObject p1, pObject p2, pObject p3);
		//!evaluates the X^2(x1,x2,x3) function, used to explore the x1, x2 (mean) space and used later approximate the minimum of the X^2 function
		/*!
			\param x1 the "mean" curvature of the positron
			\param x2 the "mean" curvature of the electron
			\param p1 the particle 1 (positron) struct
			\param p2 the particle 2 (electron) struct
			\param p3 the particle 3 (photon) struct
			\return the X^2(x1,x2,x3) value with mean curvatures x1 and x2 which also define x3(x1,x2)
		*/
		double getEliminationChiSq(double x1, double x2, pObject p1, pObject p2, pObject p3);
		
		//! Creates a contour map of the curvature space explored in the Elimination ChiSq function.  Uses the actual unsmeared event to center the contour
		/*!
			\param ap1 The particle 1 (positron) struct of the unsmeared event
			\param ap2 The particle 2 (electron) struct of the unsmeared event
			\param p1 the particle 1 (positron) struct
			\param p2 the particle 2 (electron) struct
			\param p3 the particle 3 (photon) struct
			\param bins the number of bins for the contour histogram
			\param offset the +/- range to plot in from the x1,x2 actual values, should be on the order of .001
		*/
		double generateContour(pObject ap1, pObject ap2, pObject p1, pObject p2, pObject p3, int bins, double offset);

		//! Empty constructor that initializes mass values
		DalitzChiSq();	
	
};


#endif
