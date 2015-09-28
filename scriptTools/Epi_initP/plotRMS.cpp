
#include "TFile.h"
#include <fstream>
#include <vector>
#include "TGraph.h"
#include "TCanvas.h"
#include "TMultiGraph.h"

using namespace std;
int main(){

		std::vector<double> RMSmeas;
		std::vector<double> RMSmin;
		std::vector<double> RMSelim;
		std::vector<double> initP;

		

		std::ifstream f1("RMS.txt");

		
		double temp;
		while(f1>>temp){
			RMSmeas.push_back(temp);

			f1>>temp;
			RMSmin.push_back(temp);

			f1>>temp;
			RMSelim.push_back(temp);
	
			f1>>temp;
			initP.push_back(temp);

		}
		int n = RMSmeas.size();

		double measured[n], minimalist[n], elimination[n], inp[n];

		//copy to double arrays because it would be too convenient to use vectors with root
		for(int i=0; i<n; i++){
			measured[i] = RMSmeas[i];
			minimalist[i] = RMSmin[i];
			elimination[i] = RMSelim[i];
			inp[i] = initP[i];
		}

		TCanvas *c1 = new TCanvas("c1","Epi RMS vs initial momentum",200,10,800,600);

		//TMultiGraph *mg = new TMultiGraph();

		TGraph *meas = new TGraph(n,inp,measured);
		TGraph *min = new TGraph(n,inp,minimalist);
		TGraph *elim = new TGraph(n,inp,elimination);
		
		
		

		//measured rms options
		meas->SetMarkerStyle(21);
		meas->SetMarkerColor(2);
		meas->SetLineColor(4);
		meas->Draw("ACP");
		
		min->SetMarkerColor(4);
		min->SetMarkerStyle(22);
		min->SetLineColor(3);
		min->Draw("CP");

		elim->SetMarkerColor(3);
		elim->SetMarkerStyle(3);
		elim->SetLineColor(2);
		elim->Draw("CP");

		//mg->Add(meas);
		//mg->Add(min);
		//mg->Add(elim);

		//mg->Draw("AC");

		c1->Print("E_PIrms_vs_intial_momentum.pdf");
		
	
}
