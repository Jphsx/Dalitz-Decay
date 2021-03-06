#include "PhotonHelper.h"


//These mass values should not be hardcoded
PhotonHelper::PhotonHelper(vectorFactory::ParticleParameters* event){
	evtP = event;
	M=0.13497;
	m_e=0.000511;
}
double PhotonHelper::getX3Fit(double x1, double x2){
	//DalitzChiSq* calcHelp = new DalitzChiSq();
	return mathUtility::getX3constrained(x1,x2,evtP[1],evtP[2],evtP[3]);
}
//partial derivatives for getting uncertainty on the fit photon energy
double PhotonHelper::getPartialX1(double x1, double x2){
	//DalitzChiSq* calcHelp = new DalitzChiSq();
	double a = (M*M-2*m_e*m_e);
	double b = ( (1/sin(evtP[1].theta))*(1/sin(evtP[2].theta))*mathUtility::getZ12(evtP[1],evtP[2]) );
	double c = ( (1/sin(evtP[1].theta))*mathUtility::getZ13(evtP[1],evtP[3]) );
	double d = ( (1/sin(evtP[2].theta))*mathUtility::getZ23(evtP[2],evtP[3]) );

	return ( a*c*x2*x2 + b*d ) / pow(c*x2 + d*x1,2);
}
double PhotonHelper::getPartialX2(double x1, double x2){
	//DalitzChiSq* calcHelp = new DalitzChiSq();
	double a = (M*M-2*m_e*m_e);
	double b = ( (1/sin(evtP[1].theta))*(1/sin(evtP[2].theta))*mathUtility::getZ12(evtP[1],evtP[2]) );
	double c = ( (1/sin(evtP[1].theta))*mathUtility::getZ13(evtP[1],evtP[3]) );
	double d = ( (1/sin(evtP[2].theta))*mathUtility::getZ23(evtP[2],evtP[3]) );

	return ( a*d*x1*x1 + b*c ) / pow(c*x2 + d*x1,2);
}
double PhotonHelper::getX3Variance(double x1, double x2, double x1Variance, double x2Variance, double x12covariance){
	return pow(getPartialX1(x1,x2),2)*x1Variance + pow(getPartialX2(x1,x2),2)*x2Variance + 2*getPartialX1(x1,x2)*getPartialX2(x1,x2)*x12covariance;
}
