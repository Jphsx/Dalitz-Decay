
#include "TFile.h"
#include <fstream>
#include <vector>
#include "TGraph.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
using namespace std;
int main(){

		std::vector<double> RMSmeas;
		std::vector<double> RMSmin;
		//std::vector<double> RMSelim;
		std::vector<double> initP;
		std::vector<double> RMSmin2;

		//placeholder vectors for quick and easy parsing of txt file (this isnt that efficient)
		//the initP and meas values for the plot come from minimalist1 
		std::vector<double> RMSmeas2;
		std::vector<double>initP2;

		std::vector<double>RMSmeas3;
		std::vector<double>initP3;
		

		std::ifstream f1("RMSmin.txt");
		std::ifstream f2("RMSmin2.txt");
		std::ifstream f3("RMSmeas.txt");

		
		double temp;
		while(f1>>temp){
			RMSmeas.push_back(temp);

			f1>>temp;
			RMSmin.push_back(temp);

			//f1>>temp;
			//RMSelim.push_back(temp);
	
			f1>>temp;

			initP.push_back(temp);

		}
		while(f2>>temp){
			RMSmeas2.push_back(temp);
			
			f2>>temp;
			RMSmin2.push_back(temp);
			f2>>temp;
			initP2.push_back(temp);
		}
		while(f3>>temp){
			
			RMSmeas3.push_back(temp);
			f3>>temp;
			initP3.push_back(temp);
		}
		int n = RMSmeas.size();

		double measured[n], minimalist[n], elimination[n], inp[n], minimalist2[n], measured2[n];

		//copy to double arrays because it would be too convenient to use vectors with root
		for(int i=0; i<n; i++){
			measured[i] = RMSmeas[i];
			minimalist[i] = RMSmin[i];
			//elimination[i] = RMSelim[i];
			minimalist2[i] = RMSmin2[i];
			inp[i] = initP[i];
			measured2[i] = RMSmeas3[i];
		}

		TCanvas *c2 = new TCanvas("c2","Mpi RMS vs initial momentum",200,10,800,600);

		//TMultiGraph *mg = new TMultiGraph();

		TGraph *meas = new TGraph(n,inp,measured2);
		//TGraph *min = new TGraph(n,inp,minimalist);
		//TGraph *min2 = new TGraph(n,inp,minimalist2);
		//TGraph *elim = new TGraph(n,inp,elimination);
		
		
		

		//measured rms options
		meas->SetTitle("M_{#pi^{0}} RMS vs |P_{#pi^{0}}|;Initial Momentum P_{#pi^{0}} GEV;RMS (GEV)"); 
		meas->SetMarkerStyle(21);
		meas->SetMarkerColor(2);
		meas->SetLineColor(4);
		meas->Draw("ACP");

		//meas->GetXaxis()->SetTitle("GEV");
		//meas->GetYaxis()->SetTitle("RMS");
		
	//	min->SetMarkerColor(4);
	//	min->SetMarkerStyle(22);
//		min->SetLineColor(3);
//		min->Draw("CP");

		//elim->SetMarkerColor(3);
		//elim->SetMarkerStyle(3);
		//elim->SetLineColor(2);
		//elim->Draw("CP");
		
	//	min2->SetMarkerColor(5);
	//	min2->SetMarkerStyle(20);
	//	min2->SetLineColor(5);
	//	min2->Draw("CP");

	//	TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
	//	legend->AddEntry("meas", "Measured", "p");
	//	legend->AddEntry("min", "Minimalist", "p");
	//	legend->AddEntry("min2", "Minimalist2", "p");
	//	legend->Draw();

		//mg->Add(meas);
		//mg->Add(min);
		//mg->Add(elim);

		//mg->Draw("AC");

		c2->Print("M_PIrms_vs_intial_momentum.pdf");
		
	
}
