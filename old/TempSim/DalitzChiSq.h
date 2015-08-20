#ifndef _DALITZ_CHI_SQ
#define _DALITZ_CHI_SQ
#include "TLorentzVector.h" 
#include "vectorFactory.h"
#include <cmath>
#include "TFile.h"
#include "TH2.h"


typedef vectorFactory::ParticleParameters pObject;
class DalitzChiSq
{
	
	
	public:
		double M;
		double m_e;
		double x3constrained;		

		double getMinimalistChiSq(pObject p1, pObject p2, pObject p3);
		double getEliminationChiSq(double x1, double x2, pObject p1, pObject p2, pObject p3);

		double generateContour(pObject ap1, pObject ap2, pObject p1, pObject p2, pObject p3, int bins, double offset);

		DalitzChiSq();
	private:
		
		double getVariance(pObject p);
		double getZ12(pObject p1, pObject p2);
		double getZ23(pObject p2, pObject p3);
		double getZ13(pObject p1, pObject p3);
		
		double getCosTheta(TLorentzVector v1, TLorentzVector v2);
		double getX3constrained(double x1, double x2, pObject p1, pObject p2, pObject p3);

	
};


#endif
