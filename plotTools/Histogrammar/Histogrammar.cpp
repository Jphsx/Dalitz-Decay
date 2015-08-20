#include "Histogrammar.h"
Histogrammar::Histogrammar(const char* filename, int Nevs, int modeSelect, double M, double m_e, double initialp){

	ifstream f1;
	f1.open(filename);
	if(modeSelect != -1){
	plot( f1,Nevs, modeSelect, M, m_e, initialp);
	}

}
double Histogrammar::getCosTheta(TLorentzVector v1, TLorentzVector v2){
	return (v1.Px()*v2.Px() + v1.Py()*v2.Py() + v1.Pz()*v2.Pz())/(v1.P()*v2.P());
}

void Histogrammar::plot(ifstream& fs, int Nevs, int modeSelect, double M, double m_e, double initialp){
		
		int numPs;
		double temp;
		//hardcoded sizes should be avoided
		double vdat[11];
		double bound;
		TFile* f;
		//The modes are 0 for US_CM, 1 for US_LAB, 2 for SM_LAB
		if(modeSelect==0){
			 f = new TFile("../../EventOutputs/UnsmearedCMEventHistos.root","RECREATE");
			bound =0.07;
		}
		if(modeSelect==1){ 
			f = new TFile("../../EventOutputs/UnsmearedEventHistos.root", "RECREATE");
			bound= initialp;
		}
		if(modeSelect==2){
			f = new TFile("../../EventOutputs/SmearedEventHistos.root","RECREATE");
			bound=initialp;
		}
		TLorentzVector* v = new TLorentzVector[4];
		TLorentzVector vx;
	
		
		 TH1D* hEep = new TH1D("hEep", "Energy of Positron ;Energy;Frequency of Energy",100,0,bound);
		 TH1D* hEem = new TH1D("hEem", "Energy of electron ;Energy;Frequency of Energy",100,0,bound);
		 TH1D* hEgm = new TH1D("hEgm", "Energy of photon ;Energy; Frequency of Energy",100,0,bound);
		 TH1D* hPep = new TH1D("hPep", "scalar magnitude of positron momentum ;momentum;Frequency of momentum",100,0,bound);
		 TH1D* hPem = new TH1D("hPem", "scalar magnitude of electron momentum ;momentum;Frequency of momentum",100,0,bound);
		 TH1D* hPgm = new TH1D("hPgm", "scalar magnitude of photon momentum ;momentum;Frequency of momentum",100,0,bound);
		

		 TH1D* hPhiep = new TH1D("hPhiep","Phi of positron ;phi;frequency of angle",100,-M_PI,M_PI);
		 TH1D* hPhiem = new TH1D("hPhiem","Phi of electron ;phi;frequency of angle",100,-M_PI,M_PI);
		 TH1D* hPhigm = new TH1D("hPhigm","Phi of photon ;phi;frequency of angle",100,-M_PI,M_PI);
		 TH1D* hCosep = new TH1D("hCosep","cos theta of positron ;cos;frequency of cos",100,-1.0,1.0);
		 TH1D* hCosem = new TH1D("hCosem","cos theta of electron ;cos;frequency of cos",100,-1.0,1.0);
		 TH1D* hCosgm = new TH1D("hCosgm","cos theta of photon ;cos;frequency of cos",100,-1.0,1.0);

		 TH1D* hcos12 = new TH1D("hcos12","opening angle between the electron postitron ;cos;frequency",100,-1.0,1.0);
		 TH1D* hcos13 = new TH1D("hcos13","opening angle between the electron and photon ;cos;frequency",100,-1.0,1.0);
		 TH1D* hcos23 = new TH1D("hcos23","opening angle between the positron and photon ;cos;frequency",100,-1.0,1.0);

		 TH1D* hMass = new TH1D("hMass","mass from x1+x2+x3 4 vectors;mass;frequency",100,0.12,0.14);

		// TH2D* hMpiE = new TH2D("hMpiE","Sum of Three particles energy vs Mpi invariant mass ;Energy;Mass",100,0,bound,100,0,10);
	
		 TH2D* hPhiCos12 = new TH2D("hPhiCos12","phi vs cos electron+postiron 4vector ;cos;phi",100,-1.0,1.0,100, -M_PI,M_PI);
		 TH2D* hPhiCos13 = new TH2D("hPhiCos13","phi vs cos electron+photon 4vector ;cos;phi",100,-1.0,1.0,100,-M_PI,M_PI);
		 TH2D* hPhiCos23 = new TH2D("hPhiCos23","phi vs cos positron+photon 4vector ;cos;phi",100,-1.0,1.0,100, -M_PI,M_PI);
		 TH2D* hPhiCos123 = new TH2D("hPhiCos123","phi vs cos all 3 particles 3vector ;cos;phi",100,-1.0,1.0,100,-M_PI,M_PI);	
		
		
		TH1D* hEcons = new TH1D("hEcons", "Sum of  energies;energy;frequency of Energy",100,bound-1,bound+1);
		TH1D* hPcons = new TH1D("hPcons", "Sum of  scalar magnitude of momentum;momentum;frequency of momentum",100,bound-1,bound+1);
		TH1D* hMcons = new TH1D("hMcons", "Reconstructed Mass (GeV); Mass (GeV); Events per 0.001 GeV",70,0.100,0.170);
	

		
		
		for(int z=0; z<Nevs; z++){
			if(fs.is_open()){
				fs>>numPs;

	
				for(int i = 1; i<=numPs; i++){
					for(int j=0; j<11; j++){
						fs>>temp;
						vdat[j]=temp;
						//cout<<temp<<endl;
					}
		//look at the type of particle and put it in the corresponding index, then do calculations e.g. pt
					v[i].SetXYZM(vdat[6],vdat[7],vdat[8],vdat[10]);

	
				}
				
				hEep->Fill(v[1].E());
				hEem->Fill(v[2].E());
				hEgm->Fill(v[3].E());
			
				hPep->Fill(v[1].P());
				hPem->Fill(v[2].P());
				hPgm->Fill(v[3].P());

				hPhiep->Fill(v[1].Phi());
				hPhiem->Fill(v[2].Phi());
				hPhigm->Fill(v[3].Phi());

				hCosep->Fill(v[1].CosTheta());
				hCosem->Fill(v[2].CosTheta());
				hCosgm->Fill(v[3].CosTheta());
		
				hcos12->Fill(getCosTheta(v[1],v[2]));
				hcos13->Fill(getCosTheta(v[1],v[3]));
				hcos23->Fill(getCosTheta(v[2],v[3]));
	
				vx = v[1]+v[2]+v[3];
				//hMpiE->Fill(v[1].E()+v[2].E()+v[3].E(),vx.M());
				hMass->Fill(vx.M());
				hEcons->Fill(vx.E());
				hPcons->Fill(vx.P());
				hMcons->Fill(vx.M());
				vx = v[1]+v[2];
				hPhiCos12->Fill(vx.CosTheta(),vx.Phi());
				vx = v[2]+v[3];
				hPhiCos13->Fill(vx.CosTheta(),vx.Phi());
				vx = v[1]+v[3];
				hPhiCos23->Fill(vx.CosTheta(),vx.Phi());
				vx = v[1]+v[2]+v[3];
				hPhiCos123->Fill(vx.CosTheta(),vx.Phi());
	
			}
		}

	
		f->Write();
	
}

int main(int argc, char* argv[]){
	char* filename = argv[1];
	int Nevs = atoi(argv[2]);
	int modeSelect = atoi(argv[3]);
	double M = atof(argv[4]);
	double m_e = atof(argv[5]);
	double initialp = atof(argv[6]);

	Histogrammar* h = new Histogrammar( filename, Nevs, modeSelect,  M, m_e,  initialp);
	delete h;
}

