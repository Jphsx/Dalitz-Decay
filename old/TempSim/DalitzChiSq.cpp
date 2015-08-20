#include "DalitzChiSq.h"

using namespace std;



DalitzChiSq::DalitzChiSq(){
	M = 0.13497;
	m_e = 0.000511;
}
double DalitzChiSq::getCosTheta(TLorentzVector v1, TLorentzVector v2){
	return (v1.Px()*v2.Px() + v1.Py()*v2.Py() + v1.Pz()*v2.Pz())/(v1.P()*v2.P());
}
double DalitzChiSq::getZ12(pObject p1, pObject p2){
	double beta1 = p1.v.P()/p1.v.E();
	double beta2 = p2.v.P()/p2.v.E();
	return 2*((1/(beta1*beta2)) - getCosTheta(p1.v,p2.v));
}
double DalitzChiSq::getZ23(pObject p2, pObject p3){
	double beta2 = p2.v.P()/p2.v.E();
	return 2*((1/beta2) - getCosTheta(p2.v,p3.v));
}
double DalitzChiSq::getZ13(pObject p1, pObject p3){
	double beta1 = p1.v.P()/p1.v.E();
	return 2*((1/beta1) - getCosTheta(p1.v,p3.v));
}
double DalitzChiSq::getVariance(pObject p){
	if(p.pID == 11 || p.pID == -11){
		return (2e-5*2e-5) + pow((1e-3 * 1/p.pt * 1/sin(p.theta) ),2);	
	}
	if(p.pID == 22){
		return pow(p.x_m,2) * (pow(0.16/sqrt(p.x_m),2) + (0.01*0.01));
	}
}
double DalitzChiSq::getX3constrained(double x1,double x2, pObject p1, pObject p2, pObject p3){
	double topterm =   M*M - 2*m_e*m_e - (getZ12(p1,p2)/(sin(p1.theta) * sin(p2.theta) *x1*x2));
	double bottomterm = ( (getZ13(p1,p3)/(sin(p1.theta) * x1)) + (getZ23(p2,p3)/(sin(p2.theta) * x2)) );
	return topterm/bottomterm;
}
double DalitzChiSq::getMinimalistChiSq(pObject p1, pObject p2, pObject p3){
	return pow( (getX3constrained(p1.x_m,p2.x_m,p1,p2,p3) - p3.x_m) , 2)/getVariance(p3);
}
double DalitzChiSq::getEliminationChiSq(double x1, double x2, pObject p1, pObject p2, pObject p3){
	double term1 = pow(x1-(p1.x_m),2) / getVariance(p1);
	double term2 = pow(x2-(p2.x_m),2) / getVariance(p2);
	double term3 = pow( getX3constrained(x1,x2,p1,p2,p3) - p3.x_m , 2)/getVariance(p3);
	return term1+term2+term3;
}
//ap1&ap2 is true particle values where the contour will be centered and span +/- offset with resolution by n bins
double DalitzChiSq::generateContour(pObject ap1, pObject ap2, pObject p1, pObject p2, pObject p3, int bins, double offset){
	TFile *f = new TFile("../DalitzContour.root", "update");
	int nbinsx=bins;
	int nbinsy=bins;
	
	double xmin=(ap1.x_m)-offset;
	double ymin=(ap2.x_m)-offset;
	double xmax=(ap1.x_m)+offset;
	double ymax=(ap2.x_m)+offset;
	
	double chisq;
	double x1,x2;
	
	TH2D *hcont = new TH2D("hcont","#chi^{2} contours for Dalitz event; k_{1} ; k_{2}", nbinsx, xmin, xmax, nbinsy, ymin, ymax);
	for (int i=1; i <= nbinsx; i++){                  // ROOT starts at 1 (0 is underflow)         
        	x1 = hcont->GetXaxis()->GetBinCenter(i);
        	for (int j=1; j <= nbinsy; j++){         
           		x2 = hcont->GetYaxis()->GetBinCenter(j);
           		chisq = getEliminationChiSq(x1,x2,p1,p2,p3);
           		hcont->Fill(x1,x2,chisq);
		}
	}

	f->Write();
}

/*int main(){cout<<"hello"<<endl;
DalitzChiSq* test = new DalitzChiSq();
cout<<test->m_e<<endl;}*/
