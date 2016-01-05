#include "TFile.h"
#include <fstream>
#include <iomanip>
#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include "TLorentzVector.h"
#include "../../MC_Rejection/XYgenerator.h"

#include <iomanip>
using namespace std;

int main(int argc, char* argv[]){

int N = 10000;
double M = 0.13497;
double nu = 2*0.000511/M;


XYgenerator* xy = new XYgenerator(nu);
double betaMax = xy->getBeta(nu*nu);

TFile *f = new TFile("DifferentialDecay.root","RECREATE");
// x->[1, v^2] , y->[-B,B]

TH1D* hdiffDecay = new TH1D("hdiffDecay", "differential decay rate (x,y); decay rate; event per bin", 1000, 0, 9000); 
TH2D* hxy = new TH2D("hxy","y(x) vs x ; y; x",1000,-betaMax,betaMax,1000,1,nu*nu);
TH1D* hE1 = new TH1D("hE1", "E1(x); Energy GEV; event per bin", 1000, 0, M);  
TH1D* hE2 = new TH1D("hE2", "E2(x); Energy GEV; event per bin", 1000, 0, M); 
TH1D* hE3 = new TH1D("hE3", "E3(x); Energy GEV; event per bin", 1000, 0, M); 

double E1,E2,E3;
	for(int i=0; i<N; i++){
		double* xyarr = xy->generateXYpair();
		hdiffDecay->Fill(xy->getProbability(xyarr[0],xyarr[1]));
		hxy->Fill(xyarr[0],xyarr[1]);
		hE1->Fill( M/4 * (1+xyarr[0]+xyarr[1]-xyarr[0]*xyarr[1]) );
		hE2->Fill( M/4 * (1+xyarr[0]-xyarr[1]+xyarr[0]*xyarr[1]) );
		hE3->Fill( M/2 * (1-xyarr[0]) );
	}

f->Write();

}
