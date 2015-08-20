#include <iostream>
#include "TLorentzVector.h"
#include "TMath.h"
#include "TF2.h"
#include "TH2.h"
#include "TFile.h"
#include <cmath>
#include <string>
#include <fstream>

//minimizer
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"

using namespace std;
class DalitzChiSq2{
	public:
		//main function (k1,k2,E3)->(x1,x2,x3)

		double minimalistChiSq();
		double eliminationChiSq(double x1, double x2);
		

		//calculation helpers
		double getX3Var();
		double getX12Var(double pt,double theta);

		DalitzChiSq2();

		double M;
		double m_e;

		struct ParticleParameters{
			TLorentzVector v;
			double theta;
			double pt;
			double x_m;
			
		
		};
		ParticleParameters p1; //i.e. particle 1-2-3..etc
		ParticleParameters p2;
		ParticleParameters p3;
		ParticleParameters ap1; //actual unsmeared particle 1-2-3 
		ParticleParameters ap2;
		ParticleParameters ap3;

		//gives lorentz vector to struct
		void setp1(TLorentzVector vec);
		void setp2(TLorentzVector vec);
		void setp3(TLorentzVector vec);
		void setap1(TLorentzVector vec);
		void setap2(TLorentzVector vec);
		void setap3(TLorentzVector vec);

		//makes calculations for struct based on vector
		//void evaluateParameters();

		double evaluateX3(double x1, double x2, double a12, double a23, double a13);

		//temporary stuff
		double z1,z2,z3;

		
		
};
DalitzChiSq2::DalitzChiSq2(){
	M = 0.13497;
	m_e = 0.000511;
}
double DalitzChiSq2::minimalistChiSq(){
	//em & ep
    	    double cos12 = (p1.v.Px()*p2.v.Px() + p1.v.Py()*p2.v.Py() + p1.v.Pz()*p2.v.Pz())/(p1.v.P()*p2.v.P());
    	   //em & gm
   	    double cos23 = (p2.v.Px()*p3.v.Px() + p2.v.Py()*p3.v.Py() + p2.v.Pz()*p3.v.Pz())/(p3.v.P()*p2.v.P());
            //ep & gm
            double cos13 = (p1.v.Px()*p3.v.Px() + p1.v.Py()*p3.v.Py() + p1.v.Pz()*p3.v.Pz())/(p1.v.P()*p3.v.P());

	    double z12 = 2*(1-cos12);
	    double z23 = 2*(1-cos23);
	    double z13 = 2*(1-cos13);
		

	return pow( ( evaluateX3(1/p1.pt,1/p2.pt,z12,z23,z13)-p3.v.E() ),2)/getX3Var();


}
double DalitzChiSq2::eliminationChiSq(double x1, double x2){//soon to be elimination (old minimalist)
	//em & ep
    	    double cos12 = (p1.v.Px()*p2.v.Px() + p1.v.Py()*p2.v.Py() + p1.v.Pz()*p2.v.Pz())/(p1.v.P()*p2.v.P());
    	   //em & gm
   	    double cos23 = (p2.v.Px()*p3.v.Px() + p2.v.Py()*p3.v.Py() + p2.v.Pz()*p3.v.Pz())/(p3.v.P()*p2.v.P());
            //ep & gm
            double cos13 = (p1.v.Px()*p3.v.Px() + p1.v.Py()*p3.v.Py() + p1.v.Pz()*p3.v.Pz())/(p1.v.P()*p3.v.P());
		
		double beta1 = p1.v.P()/p1.v.E();
		double beta2 = p2.v.P()/p2.v.E();
	    

	    double z12 = 2*((1/(beta1*beta2)) -cos12);
	    double z23 = 2*((1/beta2)-cos23);
	    double z13 = 2*((1/beta1)-cos13);

		z1=z12;
		z2=z23;
		z3=z13;

		//cout<<" z12 z23 z13 "<<z12<<" "<<z23<<" "<<z13<<endl;

	double term1 = pow(x1-(1/p1.pt),2) / getX12Var(p1.pt,p1.theta); //only squaring numerators since variance is sigma^2
	double term2 = pow(x2-(1/p2.pt),2) / getX12Var(p2.pt,p2.theta);
	double term3 = pow(evaluateX3(x1,x2,z12,z23,z13)-p3.v.E(),2) / getX3Var();

	return term1+term2+term3;


}
void DalitzChiSq2::setp1(TLorentzVector vec){
	p1.v.SetXYZM(vec.Px(),vec.Py(),vec.Pz(),vec.M());
	
	double Pt = TMath::Sqrt(vec.Px()*vec.Px() + vec.Py()*vec.Py());
	double p = TMath::Sqrt(Pt*Pt + vec.Pz()*vec.Pz());
	double thet = TMath::ACos(vec.Pz()/p);
	
	p1.pt=Pt;
	p1.theta=thet;
	p1.x_m=1/Pt;
	//cout<<"positron : pt theta x1m : "<<Pt<<" "<<thet<<" "<<1/Pt<<endl;
	
}
void DalitzChiSq2::setp2(TLorentzVector vec){
	p2.v.SetXYZM(vec.Px(),vec.Py(),vec.Pz(),vec.M());

	double Pt = TMath::Sqrt(vec.Px()*vec.Px() + vec.Py()*vec.Py());
	double p = TMath::Sqrt(Pt*Pt + vec.Pz()*vec.Pz());
	double thet = TMath::ACos(vec.Pz()/p);
	
	p2.pt=Pt;
	p2.theta=thet;
	p2.x_m=1/Pt;
	//cout<<"electron : pt theta x1m : "<<Pt<<" "<<thet<<" "<<1/Pt<<endl;
	
}
void DalitzChiSq2::setp3(TLorentzVector vec){
	p3.v.SetXYZM(vec.Px(),vec.Py(),vec.Pz(),vec.M());
	//cout<<"photon energy "<<vec.E()<<endl;
	
}
void DalitzChiSq2::setap1(TLorentzVector vec){
	ap1.v.SetXYZM(vec.Px(),vec.Py(),vec.Pz(),vec.M());
	
	double Pt = TMath::Sqrt(vec.Px()*vec.Px() + vec.Py()*vec.Py());
	double p = TMath::Sqrt(Pt*Pt + vec.Pz()*vec.Pz());
	double thet = TMath::ACos(vec.Pz()/p);
	
	ap1.pt=Pt;
	ap1.theta=thet;
	ap1.x_m=1/Pt;
	
}
//these can be done better, pass in a specific struct and populate that struct
//could reduce these 6 calls into 1 or 2
void DalitzChiSq2::setap2(TLorentzVector vec){
	ap2.v.SetXYZM(vec.Px(),vec.Py(),vec.Pz(),vec.M());

	double Pt = TMath::Sqrt(vec.Px()*vec.Px() + vec.Py()*vec.Py());
	double p = TMath::Sqrt(Pt*Pt + vec.Pz()*vec.Pz());
	double thet = TMath::ACos(vec.Pz()/p);
	
	ap2.pt=Pt;
	ap2.theta=thet;
	ap2.x_m=1/Pt;
	
}
void DalitzChiSq2::setap3(TLorentzVector vec){
	ap3.v.SetXYZM(vec.Px(),vec.Py(),vec.Pz(),vec.M());
	//em & ep
    	    double cos12 = (ap1.v.Px()*ap2.v.Px() + ap1.v.Py()*ap2.v.Py() + ap1.v.Pz()*ap2.v.Pz())/(ap1.v.P()*ap2.v.P());
    	   //em & gm
   	    double cos23 = (ap2.v.Px()*ap3.v.Px() + ap2.v.Py()*ap3.v.Py() + ap2.v.Pz()*ap3.v.Pz())/(ap3.v.P()*ap2.v.P());
            //ep & gm
            double cos13 = (ap1.v.Px()*ap3.v.Px() + ap1.v.Py()*ap3.v.Py() + ap1.v.Pz()*ap3.v.Pz())/(ap1.v.P()*ap3.v.P());
		
		double beta1 = ap1.v.P()/ap1.v.E();
		double beta2 = ap2.v.P()/ap2.v.E();
	    

	    double z12 = 2*((1/(beta1*beta2)) -cos12);
	    double z23 = 2*((1/beta2)-cos23);
	    double z13 = 2*((1/beta1)-cos13);
	double topterm =   M*M - 2*m_e*m_e - ( z12/(TMath::Sin(ap1.theta) * TMath::Sin(ap2.theta) *ap1.x_m*ap2.x_m));
	double bottomterm = ( (z13/(TMath::Sin(ap1.theta) * ap1.x_m)) + (z23/(TMath::Sin(ap2.theta) * ap2.x_m)) );

	ap3.x_m=topterm/bottomterm;
	
}

double DalitzChiSq2::evaluateX3(double x1, double x2, double z12, double z23, double z13){
	double topterm =   M*M - 2*m_e*m_e - ( z12/(TMath::Sin(p1.theta) * TMath::Sin(p2.theta) *x1*x2));
	double bottomterm = ( (z13/(TMath::Sin(p1.theta) * x1)) + (z23/(TMath::Sin(p2.theta) * x2)) );
	return topterm/bottomterm;
}
double DalitzChiSq2::getX3Var(){
	return (p3.v.E()*p3.v.E())*( (0.16/TMath::Sqrt(p3.v.E()))*(0.16/TMath::Sqrt(p3.v.E())) + (0.01*0.01) );
}
double DalitzChiSq2::getX12Var(double pt,double theta){
	return  (2e-5*2e-5) + pow((1e-3 * 1/pt * 1/TMath::Sin(theta) ),2);
}



int main(){
//the unsmeared sample 4vectors
//1 22 0 0 0 0 -4.1078247243 1.68355857009 0.790457556728 4.50925898274 5.96046447754e-08
//1 11 0 0 0 0 -2.57783720676 1.1352928704 0.475596305383 2.85662855005 0.000510999998632
//1 -11 0 0 0 0 -2.37767887022 1.04769700879 0.438544803129 2.63502327079 0.000511000002977

//sample 4 vector, this one is already smeared should have unsmeared and call a smearing class
//1 22 0 0 0 0 -4.71429479708 1.93211539932 0.907158945939 5.17499591336 0.0 -0.0 0.0 0.0 0.0
//1 11 0 0 0 0 -2.57276133106 1.13305742843 0.474659835181 2.85100372202 0.000510999998632 -0.0 0.0 0.0 0.0
//1 -11 0 0 0 0 -2.37299823439 1.04563454013 0.437681495412 2.62983603377 0.000511000002977 -0.0 0.0 0.0 0.0

	
	TFile *f = new TFile("Dalitz_chi.root", "update");

	TLorentzVector v3;
	TLorentzVector v2;
	TLorentzVector v1;
	//smeared
	v3.SetXYZM(-4.71429479708, 1.93211539932, 0.907158945939, 0.0);
	v2.SetXYZM(-2.57276133106, 1.13305742843, 0.474659835181, 0.000510999998632);
	v1.SetXYZM(-2.37299823439, 1.04563454013, 0.437681495412, 0.000511000002977);

	//unsmeared
	TLorentzVector av3;
	TLorentzVector av2;
	TLorentzVector av1;
	
	av3.SetXYZM(-4.1078247243, 1.68355857009, 0.790457556728, 5.96046447754e-08);
	av2.SetXYZM(-2.57783720676, 1.1352928704, 0.475596305383, 0.000510999998632);
	av1.SetXYZM(-2.37767887022, 1.04769700879, 0.438544803129, 0.000511000002977);

	DalitzChiSq2 *dcs = new DalitzChiSq2();
	dcs->setp1(v1);
	dcs->setp2(v2);
	dcs->setp3(v3);
	dcs->setap1(av1);
	dcs->setap2(av2);
	dcs->setap3(av3);

	
	
	cout<<"The Real x1, x2, x3 values: "<<1/dcs->ap1.pt<<" "<<1/dcs->ap2.pt<<" "<<dcs->ap3.x_m <<endl;
	cout<<"The Measured Values x1m x2m x3m: "<<1/dcs->p1.pt<<" "<<1/dcs->p2.pt<<" "<<dcs->p3.v.E()<<endl;
	//cout<<"The fit values x1f, x2f, x3f: "<<0.385634<<" "<<0.355722<<" "<<dcs->evaluateX3(0.385634,0.355722,dcs->z1,dcs->z2,dcs->z3)<<endl;
	//cout<<"zs "<<dcs->z1<<" "<<dcs->z2<<" "<<dcs->z3<<endl;
	
	//x         = 0.385634     +/-  0.000391586
	//y         = 0.355722     +/-  0.000361297

	//starting with hardcoded bins
	int nbinsx=200;
	int nbinsy=200;
	double xmin=(1/dcs->ap1.pt)-.001;
	double ymin=(1/dcs->ap2.pt)-.001;
	double xmax=(1/dcs->ap1.pt)+.001;
	double ymax=(1/dcs->ap2.pt)+.001;

	double chisq;
	double x1,x2;

	   TH2D *h2 = new TH2D("h2","#chi^{2} contours for Dalitz event; k_{1} ; k_{2}", nbinsx, xmin, xmax, nbinsy, ymin, ymax);

   for (int i=1; i <= nbinsx; i++){                  // ROOT starts at 1 (0 is underflow)         
        x1 = h2->GetXaxis()->GetBinCenter(i);
        for (int j=1; j <= nbinsy; j++){         
           x2 = h2->GetYaxis()->GetBinCenter(j);
           chisq = dcs->eliminationChiSq(x1,x2); //call function
           h2->Fill(x1,x2,chisq);
//           outfile << vals[0] << " " << vals[1] << " " << chisqb << endl; 
        }
   }
	//cout<<"CHISQ: "<<dcs->minimalistChiSq();
	cout<<"The fit values x1f, x2f, x3f: "<<0.385634<<" "<<0.355722<<" "<<dcs->evaluateX3(0.385634,0.355722,dcs->z1,dcs->z2,dcs->z3)<<endl;
   //h2->Draw();
	f->Write();


}
