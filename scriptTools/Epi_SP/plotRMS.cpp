
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
		std::vector<double> RMSelim;
		std::vector<double> scaleParameterNP;

		

		std::ifstream f1("RMS.txt");

		
		double temp;
		while(f1>>temp){
			RMSmeas.push_back(temp);

			f1>>temp;
			RMSmin.push_back(temp);

			f1>>temp;
			RMSelim.push_back(temp);
	
			f1>>temp;
			scaleParameterNP.push_back(temp);

		}
		int n = RMSmeas.size();

		double measured[n], minimalist[n], elimination[n], sp[n];

		//copy to double arrays because it would be too convenient to use vectors with root
		for(int i=0; i<n; i++){
			measured[i] = RMSmeas[i];
			minimalist[i] = RMSmin[i];
			elimination[i] = RMSelim[i];
			sp[i] = scaleParameterNP[i];
		}

		TCanvas *c1 = new TCanvas("c1","Epi RMS vs scaleparameter",200,10,800,600);

		TMultiGraph *mg = new TMultiGraph();

		TGraph *meas = new TGraph(n,sp,measured);
		TGraph *min = new TGraph(n,sp,minimalist);
		TGraph *elim = new TGraph(n,sp,elimination);
		
		mg->SetTitle("E_#pi^{0} RMS vs #sigma_sp;#sigma;RMS"); 
		

		//measured rms options
		meas->SetMarkerStyle(21);
		meas->SetMarkerColor(2);
		meas->SetLineColor(4);
		meas->SetTitle("Measured");
		//meas->Draw("ACP");
		
		min->SetMarkerColor(4);
		min->SetMarkerStyle(22);
		min->SetLineColor(3);
		min->SetTitle("Minimalist");
		//min->Draw("CP");

		elim->SetMarkerColor(3);
		elim->SetMarkerStyle(3);
		elim->SetLineColor(2);
		elim->SetTitle("Elimination");
		//elim->Draw("CP");

		mg->Add(meas);
		mg->Add(min);
		mg->Add(elim);

		
		/*TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
		legend->AddEntry("meas", "Measured");
		legend->AddEntry("min", "Minimalist");
		legend->AddEntry("elim", "Elimination");
		//legend->Draw();
		mg->Add(legend);*/

		mg->Draw("ACP");
		c1->BuildLegend();

		c1->Print("E_PIrms_vs_sigma_sp.pdf");
		
	
}
