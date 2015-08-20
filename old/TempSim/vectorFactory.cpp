#include "vectorFactory.h"
using namespace std;

void vectorFactory::readVector(const char* path){
	
	ifstream f (path);
	
	int numPs;
	double temp;
	//hardcoded sizes should be avoided
	double vdat[11];
	
	

	if(f.is_open()){
	f>>numPs;

	ParticleParameters *ptr = new ParticleParameters [numPs+1];
	arrPtr = ptr;
	
	for(int i = 0; i<numPs; i++){
		for(int j=0; j<11; j++){
			f>>temp;
			vdat[j]=temp;
			//cout<<temp<<endl;
		}
		//look at the type of particle and put it in the corresponding index, then do calculations e.g. pt
		
		setTLVector((int)vdat[1],vdat[6],vdat[7],vdat[8],vdat[10]);
	
	}
	
	}


}
void vectorFactory::setTLVector(int pID,double px,double py,double pz, double m){
	if(pID == -11){
		arrPtr[1].v.SetXYZM(px,py,pz,m);
		arrPtr[1].pID=pID;
		makeCalculations(arrPtr[1]);
	}
	if(pID == 11){
		arrPtr[2].v.SetXYZM(px,py,pz,m);
		arrPtr[2].pID=pID;
		makeCalculations(arrPtr[2]);
	}
	if(pID == 22){
		arrPtr[3].v.SetXYZM(px,py,pz,m);
		arrPtr[3].pID=pID;
		makeCalculations(arrPtr[3]);
	}
	

}
void vectorFactory::makeCalculations(ParticleParameters& particle){
	if(particle.pID==11 || particle.pID==-11){
		double Pt = TMath::Sqrt(particle.v.Px()*particle.v.Px() + particle.v.Py()*particle.v.Py());
		double p = TMath::Sqrt(Pt*Pt + particle.v.Pz()*particle.v.Pz());
		double thet = TMath::ACos(particle.v.Pz()/p);
	
		particle.pt=Pt;
		particle.theta=thet;
		particle.x_m=1/Pt;
	}
	//else return;
	if(particle.pID==22){
		particle.pt=0.0;
		particle.theta=0.0;
		particle.x_m=particle.v.E();
	}
	
}
/*int main(){
	vectorFactory *v = new vectorFactory();

	v->readVector("dalitz_0mm_unsmeared.hepevt");
	cout<<"about to print"<<endl;
	v->printVectors();
	
}*/
