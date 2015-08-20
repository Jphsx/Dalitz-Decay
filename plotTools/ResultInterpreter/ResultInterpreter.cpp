
#include "ResultInterpreter.h"
ResultInterpreter::ResultInterpreter(const char* file1, int Nevs){
	ifstream f1;
	f1.open(file1);
	makeHistos(Nevs, f1);

}
void ResultInterpreter::makeHistos(int Nevs, ifstream& f1){
	int numPs;
	double temp;
	//hardcoded sizes should be avoided
	const int SIZE = 6;
	double fitdat[SIZE];
	//fit dat { xm , xm stdev, x actual, xfit, xfit stdev}
	double pullVal;
	double pullValrf;
	TFile *f = new TFile("../../EventOutputs/ResultHistos.root", "UPDATE");

	//initialize pull distributions (has to be done localy TH1 cannot be accessed as a class member)
	TH1D* hpullx1 = new TH1D("hpullx1", "x1 pull distribution;(xfit-xm)/sqrt(vm-vfit);Frequency",100,-5,5);
	TH1D* hpullx2 = new TH1D("hpullx2", "x2 pull distribution;(xfit-xm)/sqrt(vm-vfit);Frequency",100,-5,5);
	TH1D* hpullx3 = new TH1D("hpullx3", "x3 pull distribution;(xfit-xm)/sqrt(vm-vfit);Frequency",100,-5,5);
	TH1D* hpullrfx1 = new TH1D("hpullrfx1", "pull distribution between x1 real and x1 fit; x1_real-x1fit/varx2;Frequency",100,.1,0.3);
	TH1D* hpullrfx2 = new TH1D("hpullrfx2", "pull distribution between x2 real and x2 fit; x2_real-x2fit/varx2;Frequency",100,-0.2,.1);
	TH1D* hpullrfx3 = new TH1D("hpullrfx3", "pull distribution between x3 real and x3 fit; x3_real-x3fit/varx3;Frequency",100,-0.3,0.3);
	TH1D* hchisq = new TH1D("hchisq", " elimination #chi^{2}; #chi^{2}; Frequency",100,0,10);
	TH1D* hchisqLarge = new TH1D("hchisqLarge","elimination #chi^{2};#chi^{2};Events Per bin",100,0,30);
	



	if(f1.is_open()){
		for(int z=0; z<Nevs; z++){
			f1>>numPs;
			//cout<<numPs<<endl;
	
			for(int i = 0; i<numPs; i++){
				for(int j=0; j<SIZE; j++){
				
					f1>>temp;
					
					//cout<<temp<<endl;
					fitdat[j]=temp;
						
					
				}
				//fill pull distributions
				pullVal = (fitdat[3]-fitdat[0])/sqrt(fitdat[1]*fitdat[1] - fitdat[4]*fitdat[4]);
				pullValrf = (fitdat[2]-fitdat[3])/fitdat[4];
				//cout<< pullVal<<endl;
				if(i==0) {
					hpullx1->Fill(pullVal);
					hpullrfx1->Fill(pullValrf);
				}
				if(i==1){
					 hpullx2->Fill(pullVal);
					 hpullrfx2->Fill(pullValrf);
				}
				if(i==2) {
					 hpullx3->Fill(pullVal);
					 hpullrfx3->Fill(pullValrf);
				}
				if(i==0) {
					hchisq->Fill(fitdat[5]);
					hchisqLarge->Fill(fitdat[5]);
				}

				
			}
		}
	}
	f->Write();
	
}
int main( int argc, char* argv[] ){
	const char* filename = argv[1];
	int Nevs = atoi(argv[2]);
	ResultInterpreter* r = new ResultInterpreter(filename, Nevs);
	delete r;
}
