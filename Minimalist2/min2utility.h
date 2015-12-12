#ifndef _MIN2HELP_
#define _MIN2HELP_
#include "TLorentzVector.h"
#include <cmath>
#include "../mathUtility/mathUtility.h"
#include <iostream>
#include <iomanip>
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMultiGraph.h"

using namespace std;
class min2utility{
	
	public:
	//method to compute X^2
	//double getChiSqMin2(double Mpi, double M12, double E12, double P12, double psi12_3);

	//overloaded X^2 that takes TLorentz Vector
	double getChiSqMin2(double Mpi, TLorentzVector v12, TLorentzVector v3, double psi123, double psi123m);

	//bisection to minimize derivatives with given measured values
	
	//method for derivatives

	//method for computing variance based on given scale parameter for sigma_psi_12,3 ( this is actually trivial and in constructor)
	
	min2utility(double M, TLorentzVector v12, TLorentzVector v3);
	min2utility(double M);
	//sets the vectors so only 1 instance is needed in minimalist
	void setVectors(TLorentzVector v12, TLorentzVector v3);
	//static double bisection(double(min2utility::*f)(double), double a, double b, double TOL, int N );
	double bisection(double a, double b, double TOL, int N);
	double sign(double param);
	//NOTE: use mathutility to get x3constrainedMin2
	
	//variables to store locally v12, and v3?	
	//double scaleParameterCP;
	double var_psi12_3;
	TLorentzVector x12_u;
	TLorentzVector x3_u;
	double Mpi;
	
	double MinimizeMin2();
	
	//helper derivative method for minimization
	double Min2ChiPrime(double psi123);

};
#endif
