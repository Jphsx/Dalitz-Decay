
//Here is the actual event to be copied, Event #3 from 10GEV initP
//1 -11 0 0 0 0 -0.951264672 1.56248772 -0.283606968 1.85113628 0.000511
//1 11 0 0 0 0 -0.224975998 0.370413795 -0.0654394757 0.438295751 0.000511
//1 22 0 0 0 0 -3.75269207 6.61996006 -1.24913409 7.71147878 1.1920929e-07

#include "TFile.h"
#include <fstream>
#include <iomanip>
#include "TH1D.h"
#include <iostream>
#include "TLorentzVector.h"
#include "../../Smearing/DalitzSmear.h"
#include "../../mathUtility/mathUtility.h"
#include <iomanip>
using namespace std;

int main(int argc, char* argv[]){

	//number of iterations to run
	//int N = argv[1];
	int N = 10000;
	
TLorentzVector v12, v1,v2,v3, v12sm, v1sm, v2sm, v3sm;

v1.SetXYZM(-0.951264672, 1.56248772, -0.283606968, 0.000511);
v2.SetXYZM(-0.224975998, 0.370413795, -0.0654394757, 0.000511);
v3.SetXYZM( -3.75269207, 6.61996006, -1.24913409, 1.1920929e-07);
 
v12 = v1+v2;
double spCP = 1e-6;
double spNP = 1e-4;

TFile *f = new TFile("OpeningAngleVariance.root","RECREATE");
TH1D* hpsi123 = new TH1D("hpsi123", "distribution of opening angles of a 10GEV event; psi123(radians); event per bin",100,0,0.25);
TH1D* hpsiDiff = new TH1D("hpsiDiff", "Actual and Smeared opening angle Difference; (rads); event per bin", 100, -0.25, 0.25);

//checks to see if both directions are normally distributed
TH1D* hAngSm1theta = new TH1D("hAngSm1theta", "particle1 angular smearing on theta; (rads); event per bin", 100, 1.72460907-2*spCP, 1.72460907+2*spCP);
TH1D* hAngSm2theta = new TH1D("hAngSm2theta", "particle2 angular smearing on theta; (rads); event per bin", 100, 1.72066116-2*spCP, 1.72066116+2*spCP);
TH1D* hAngSm3theta = new TH1D("hAngSm3theta", "particle3 angular smearing on theta; (rads); event per bin", 100, 1.73349693-2*spNP, 1.73349693+2*spNP);
TH1D* hAngSm12theta = new TH1D("hAngSm12theta", "particle12 angular smearing on theta; (rads); event per bin", 100, 1.72385329-2*spCP, 1.72385329+2*spCP);

TH1D* hAngSm1phi = new TH1D("hAngSm1phi", "particle1 angular smearing on phi; (rads); event per bin", 100,  2.11767164-2*spCP, 2.11767164+2*spCP);
TH1D* hAngSm2phi = new TH1D("hAngSm2phi", "particle2 angular smearing on phi; (rads); event per bin", 100, 2.11661293-2*spCP, 2.11661293+2*spCP);
TH1D* hAngSm3phi = new TH1D("hAngSm3phi", "particle3 angular smearing on phi; (rads); event per bin", 100, 2.08650327-2*spNP, 2.08650327+2*spNP);
TH1D* hAngSm12phi = new TH1D("hAngSm12phi", "particle12 angular smearing on phi; (rads); event per bin", 100, 2.11746886-2*spCP, 2.11746886+2*spCP);

//check to see if dpsi is actually rayleigh distributed

TH1D* hDpsi1 = new TH1D("hDpsi1", "particle1 dpsi; (rads); event per bin", 100, 0, 1e-5);
TH1D* hDpsi2 = new TH1D("hDpsi2", "particle1 dpsi; (rads); event per bin", 100, 0, 1e-5);
TH1D* hDpsi3 = new TH1D("hDpsi3", "particle1 dpsi; (rads); event per bin", 100, 0, 1e-3);

double openingActual, openingSmear;
//dpsi checks
double opening1,opening2,opening3;
for(int i=0; i<N; i++){

     	 DalitzSmear* sm = new DalitzSmear();
        //smear 1
        sm->setpid(-11);
        v1sm = sm->smearDirection(v1);
        //smear 2
        sm->setpid(11);
        v2sm = sm->smearDirection(v2);
        //smear 3
        sm->setpid(22);
        v3sm = sm->smearDirection(v3);

	v12sm = v1sm + v2sm;

	//checks framework
	hAngSm1theta->Fill(v1sm.Theta());
	hAngSm2theta->Fill(v2sm.Theta());
	hAngSm3theta->Fill(v3sm.Theta());
	hAngSm12theta->Fill(v12sm.Theta());

	
	hAngSm1phi->Fill(v1sm.Phi());
	hAngSm2phi->Fill(v2sm.Phi());
	hAngSm3phi->Fill(v3sm.Phi());
	hAngSm12phi->Fill(v12sm.Phi());

	opening1 = mathUtility::safeAcos(mathUtility::getCosTheta(v1,v1sm));
	opening2 = mathUtility::safeAcos(mathUtility::getCosTheta(v2,v2sm));
	opening3 = mathUtility::safeAcos(mathUtility::getCosTheta(v3,v3sm));

	hDpsi1->Fill(opening1);
	hDpsi2->Fill(opening2);
	hDpsi3->Fill(opening3);
	//end checks

	
	
	openingActual = mathUtility::safeAcos(mathUtility::getCosTheta(v12,  v3));
	openingSmear = mathUtility::safeAcos(mathUtility::getCosTheta( v12sm, v3sm));


	hpsi123->Fill(openingSmear);
	hpsiDiff->Fill(openingActual-openingSmear);

}
f->Write();
cout<<setprecision(9);
cout<<"True Values:"<<endl;
cout<<"p1(r,theta,phi) -> theta: "<<v1.Theta()<<" phi: "<<v1.Phi()<<endl;
cout<<"p2(r,theta,phi) -> theta: "<<v2.Theta()<<" phi: "<<v2.Phi()<<endl;
cout<<"p3(r,theta,phi) -> theta: "<<v3.Theta()<<" phi: "<<v3.Phi()<<endl;
cout<<"p12(r,theta,phi) -> theta: "<<v12.Theta()<<" phi: "<<v12.Phi()<<endl;
cout<<"Opening Angle RMS: "<< hpsi123->GetRMS()<<endl;
}



