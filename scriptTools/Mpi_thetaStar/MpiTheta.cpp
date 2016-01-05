#include "../../Minimization/MinHelper.h"
//#include "../../Minimization/PhotonHelper.h"
#include <iostream>
#include <cmath>
#include "TLorentzVector.h"
#include <fstream>
#include "TFile.h"
#include "TH2D.h"
#include "TProfile.h"
#include "../../mathUtility/mathUtility.h"
#include <iomanip>


using namespace std;

int main(int argc, char* argv[]){
int N  = atoi(argv[1]);
double M = atof(argv[2]);
double m_e = atof(argv[3]);
double initP = atof(argv[4]);
TFile *f = new TFile("dmThetaStar.root", "RECREATE");
TH2D* hdmTheta = new TH2D("hdmTheta", "Relative Mass Uncertainty vs detector #pi^{0} and detector X_{3} angle;cos#theta*; #frac{#Delta M}{M} ",1000,-1,1 ,1000,-1 ,1 );
//TH2D* hdmGenTheta = new TH2D("hdmGenTheta","Relative Mass Uncertainty vs generator #pi^{0} and detector X_{3} angle; #theta *; #frac{#Delta M}{M}",1000,0,3.14,1000,0,1);
TProfile* havgdm = new TProfile("havgdm", "Expected Relative Mass Uncertainty vs detector #pi^{0} and CM X_{3} angle; cos#theta*; #frac{#Delta M}{M}",20, -1,1,-1,1,"");

TProfile* hrmsdm = new TProfile("hrmsdm", "Relative Mass Uncertainty RMS vs detector #pi^{0} and CM X_{3} angle; cos#theta*;#frac{#Delta M}{M}",20,-1,1,-1,1,"S");

TH2D* hE3labthetastar = new TH2D("hE3labthetastar", "unsmeared photon lab energy vs lab #pi^{0} and CM X_{3} angle;cos#theta*; E3lab ",50,-1,1 ,100,0 ,10 );

MinHelper* h = new MinHelper("../../EventOutputs/DalitzSmearVectors.hepevt","../../EventOutputs/DalitzCMVectors.hepevt");
MinHelper* h2 = new MinHelper("../../EventOutputs/DalitzActualVectors.hepevt","../../EventOutputs/DalitzCMVectors.hepevt");

TLorentzVector labActual4vector;
TLorentzVector Mpi;
double costhetaStar, dM, costhetaStarMgen;

//actual vectors -> evtP , CM vectors -> evtP_actual in h2
//for a single ensemble
for(int i=0; i<N; i++){

		h->populateParticles(h->f1,h->factory1,h->evtP);
		h->populateParticles(h->f2,h->factory2,h->evtP_actual);

		h2->populateParticles(h2->f1,h2->factory1,h2->evtP);
		h2->populateParticles(h2->f2,h2->factory2,h2->evtP_actual);

		labActual4vector = h2->evtP[1].v + h2->evtP[2].v + h2->evtP[3].v ;
		Mpi = h->evtP[1].v + h->evtP[2].v + h->evtP[3].v ;

		costhetaStarMgen = mathUtility::getCosTheta(labActual4vector, h2->evtP_actual[3].v);
		costhetaStar = mathUtility::getCosTheta(Mpi,h->evtP_actual[3].v) ;
		dM = M - Mpi.M();

		//thetaStarMgen = mathUtility::safeAcos(mathUtility::getCosTheta(Mpi_actual,h->evtP[3].v) );
		hE3labthetastar->Fill(costhetaStarMgen , h2->evtP[3].v.E() );
		hdmTheta->Fill(costhetaStarMgen, dM/M);
		havgdm->Fill(costhetaStarMgen, dM/M);
		hrmsdm->Fill(costhetaStarMgen, dM/M);
		//hdmGenTheta->Fill(dM/M, thetaStarMgen);

		

}
f->Write();


}
