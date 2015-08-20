#ifndef _INTERPRETER_
#define _INTERPRETER_
#include "TFile.h"
#include "TH1D.h"
#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;


class ResultInterpreter{
	public:
		
		
		ResultInterpreter(const char* file1, int Nevs);
		void makeHistos(int Nevs, ifstream& f1);
		
	
};
#endif