#include <iostream>
#include "TLorentzVector.h"
#include "TMath.h"
#include <cmath>
#include <string>
#include <fstream>
using namespace std;

class DalitzChiSq{
	public:
		//initialization and file i/o
		//ifstream startFileStream(string filename);
		bool readVectors(ifstream &f);//maybe return pointer to array of vector values

		//math functions
		double evaluateE3(double M, double m);
		double evaluateChiSqMini();
		double evaluateChiSqElim(double k1, double k2, double k1m, double k2m, double vark1, double vark2, double varE3);
		double calculateVarCP(TLorentzVector v);
		double calculateVarNP(TLorentzVector v);
		
		//plotting and output
		//void initializePlots();
		//void populatePlots();
		//void writeOutput();
		
		//constructor and initialization
		DalitzChiSq(int events);
		//void getMeanCurvatures(string filename);

			//data structures that will house input from file
			double vdat1[15];//positron
			double vdat2[15];//electron
			double vdat3[15];//photon
		
			//vectors that encapsulate important parameters from data input
			TLorentzVector v1;//positron
			TLorentzVector v2;//electron
			TLorentzVector v3;//photon
                
		//get measured curvature
		double getKm(TLorentzVector v);
		
		//average curvatures
		double k1;
		double k2;
		
		//Number of events, for debugging and checking purposes
		int Nevents;
		void printNevents();
                
		//particle masses
		double m_pion;                   
		double m_e;
	
};
//constructor, initializes all the histograms for filling
DalitzChiSq::DalitzChiSq(int events){
	Nevents=events;
	m_pion = 0.13497;
	m_e = 0.000511;
	

}
//first function that should be called, opens up the 4vector file and grabs the k from a .root file
/*ifstream DalitzChiSq::startFileStream(string filename,string rootfile){ // plan on just passing fstream to main because why not
	getMeanCurvatures(rootfile);
	ifstream f;
	f.open(filename);
	return &f;
}*/
//reads all the vector data into arrays and then populates the corresponding TLorentzVectors with the data
bool DalitzChiSq::readVectors(ifstream &f){
	double temp;
	f>>temp;
	
	if(temp==3 && !f.eof()){
		for(int i =0; i<15;i++){
			f>>temp;
			vdat3[i]=temp;
		}
		v3.SetXYZM(vdat3[6],vdat3[7],vdat3[8],vdat3[10]);
		for(int i =0; i<15;i++){
			f>>temp;
			vdat2[i]=temp;
		}
		v2.SetXYZM(vdat2[6],vdat2[7],vdat2[8],vdat2[10]);
		for(int i= 0; i<15;i++){
			f>>temp;
			vdat1[i]=temp;
		}
		v1.SetXYZM(vdat1[6],vdat1[7],vdat1[8],vdat1[10]);
		Nevents++;
		return true;
		
	}
	else return false;	
	
	
}
//opens up a .root file from the argument, takes the mean from those histograms
/*void DalitzChiSq::getMeanCurvatures(string filename){//more specifically it takes kappa value from the kappa curve, generated by smeared 4vectors for electron and positron (uses histos_dalitz.root)
	
	TFile f(rootfile);
	TH1D *hk1 = (TH1D*)f.Get("hkapem");
	TH1D *hk2 = (TH1D*)f.Get("hkapep");

        //Bool_t flag = kTRUE;
	//hk1->StatOverflows(flag);
	//hk2->StatOverflows(flag);	

	k1 = hk1->getMean();
	k2 = hk2->getMean();
	
}*/
//returns the photon energy recalculated after mass constraint using Eq 7.3 in documentation
double DalitzChiSq::evaluateE3(double M, double m){ 
	//only exectuted after a call to readVectors
	
	    //em & ep
    	    double cos12 = (v1.Px()*v2.Px() + v1.Py()*v2.Py() + v1.Pz()*v2.Pz())/(v1.P()*v2.P());
    	   //em & gm
   	    double cos23 = (v2.Px()*v3.Px() + v2.Py()*v3.Py() + v2.Pz()*v3.Pz())/(v3.P()*v2.P());
            //ep & gm
            double cos13 = (v1.Px()*v3.Px() + v1.Py()*v3.Py() + v1.Pz()*v3.Pz())/(v1.P()*v3.P());

    		double z12 = 2.0*(1.0-cos12);
    		double z13 = 2.0*(1.0-cos13);
    		double z23 = 2.0*(1.0-cos23);
	
	 return (M*M - 2*m*m - v1.P()*v2.P()*z12) / ( (v1.P()*z13) + (v2.P()*z23) );
}
//evalutes the X^2 for an event using the minimalist method, i.e. evaluates Eq 7.5 in documentation
double DalitzChiSq::evaluateChiSqMini(){
	return pow( (evaluateE3(m_pion,m_e)-v3.E())/calculateVarNP(v3),2 );
	
}
double DalitzChiSq:: evaluateChiSqElim(double k1, double k2, double k1m, double k2m, double vark1, double vark2, double varE3){
	return 1.0;
}
//returns the variance of Charged particle, in this case either the electron or positron variance
double DalitzChiSq::calculateVarCP(TLorentzVector v){

	double Pt = TMath::Sqrt(v.Px()*v.Px() + v.Py()*v.Py());
	double p = TMath::Sqrt(Pt*Pt + v.Pz()*v.Pz());
	double theta = TMath::ACos(v.Pz()/p);
	return  (2e-5*2e-5) + ( (1e-3 * 1/Pt * 1/TMath::Sin(theta) )*(1e-3 * 1/Pt * 1/TMath::Sin(theta) ));

}
//returns the variance of Neutral Particle, in this case the photon
double DalitzChiSq::calculateVarNP(TLorentzVector v){
	return (v.E()*v.E())*( (0.16/TMath::Sqrt(v.E()))*(0.16/TMath::Sqrt(v.E())) + (0.01*0.01) );
}
double DalitzChiSq::getKm(TLorentzVector v){
	return 1/TMath::Sqrt(v.Px()*v.Px() + v.Py()*v.Py());
}	

////////////////////////////////////////////////////////

int main(){
	TH2D *hkkContMin = new TH2D("hkkContMin","curvature weighted by X^2;k1;k2",100,0.0,10.0,100,0.0,10.0);

	DalitzChiSq *dcs = new DalitzChiSq(1);
	TFile rootf = new TFile("dalitz_chi.root", "update");
	ifstream f;//= startFileStream("dalitz0mm.hepevt","histos_dalitz.root");
	f.open("dalitz0mm.hepevt");

	while(readVectors(&f)){
		hkkContMin->Fill(dcs->getKm(dcs->v1),dcs->getKm(dcs->v2),dcs->evaluateChiSqMini());
	}



	rootf.Write();
}


