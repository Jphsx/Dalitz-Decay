#ifndef _HISTOS
#define _HISTOS

#include <math.h>
#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h" 
#include <iostream>
#include <fstream>
#include <cmath>


using namespace std;

class Histogrammar{
	public:
		void plot(ifstream& fs, int Nevs, int modeSelect, double M, double m_e, double initialp);

		//mode select constructor
		Histogrammar(const char* filename, int Nevs, int modeSelect, double M, double m_e, double initialp);
		
		double getCosTheta(TLorentzVector v1, TLorentzVector v2);


};
#endif
