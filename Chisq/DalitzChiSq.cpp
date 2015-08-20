#include "DalitzChiSq.h"

using namespace std;


//constructor that initializes the masses for most calculations, need to add another one that sets these values externally
DalitzChiSq::DalitzChiSq(){
	M = 0.13497;
	m_e = 0.000511;
	
}

//evaluates the X^2(x1,x2,x3) function such that x1 x2 are exactly 0 and the mean x3 is given by the mass constrained energy
double DalitzChiSq::getMinimalistChiSq(pObject p1, pObject p2, pObject p3){
	return pow( (mathUtility::getX3constrained(p1.x_m,p2.x_m,p1,p2,p3) - p3.x_m) , 2)/mathUtility::getVariance(p3);
}
//evaluates the X^2(x1,x2,x3) function, used to explore the x1, x2 (mean) space and later approximate the minimum of the X^2 function
double DalitzChiSq::getEliminationChiSq(double x1, double x2, pObject p1, pObject p2, pObject p3){
	double term1 = pow(x1-(p1.x_m),2) / mathUtility::getVariance(p1);
	double term2 = pow(x2-(p2.x_m),2) / mathUtility::getVariance(p2);
	double term3 = pow( mathUtility::getX3constrained(x1,x2,p1,p2,p3) - p3.x_m , 2)/mathUtility::getVariance(p3);
	return term1+term2+term3;
}
//ap1&ap2 is true particle values where the contour will be centered and span +/- offset with resolution by n bins
//sets up the conditions for contour mapping of the elimination Chisq method, spans the x1,x2 space and creates a contour
double DalitzChiSq::generateContour(pObject ap1, pObject ap2, pObject p1, pObject p2, pObject p3, int bins, double offset){
	TFile *f = new TFile("../EventOutputs/DalitzContour.root", "UPDATE");
	int nbinsx=bins;
	int nbinsy=bins;
	
	double xmin=(ap1.x_m)-offset;
	double ymin=(ap2.x_m)-offset;
	double xmax=(ap1.x_m)+offset;
	double ymax=(ap2.x_m)+offset;
	
	double chisq;
	double x1,x2;
	
	TH2D *hcont = new TH2D("h2","#chi^{2} contours for Dalitz event; k_{1} ; k_{2}", nbinsx, xmin, xmax, nbinsy, ymin, ymax);
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
