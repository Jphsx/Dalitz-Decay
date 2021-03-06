#include "MinHelper.h"

MinHelper::MinHelper(const char* file1, const char* file2){
	//populateParticles(f1,f2);
	f1.open(file1);
	f2.open(file2);
	factory1 = new vectorFactory();
	factory2 = new vectorFactory();
	//populateParticles(f1,factory1,evtP);
	//populateParticles(f2,factory2,evtP_actual);
}
/*void MinHelper::populateParticles(const char* f1, const char* f2){
	vectorFactory* genSmear = new vectorFactory();
	genSmear->readVector(f1);
	evtP=genSmear->arrPtr;
	vectorFactory* genActual = new vectorFactory();
	genActual->readVector(f2);
	evtP_actual=genActual->arrPtr;
}*/
void MinHelper::populateParticles(ifstream& fs,vectorFactory* vf,vectorFactory::ParticleParameters*& evtp){

	
	int numPs;
	double temp;
	//hardcoded sizes should be avoided
	double vdat[11];
	
	

	if(fs.is_open()){
		fs>>numPs;

		vectorFactory::ParticleParameters *ptr = new vectorFactory::ParticleParameters [numPs+1];
		vf->arrPtr = ptr;
	
		for(int i = 0; i<numPs; i++){
			for(int j=0; j<11; j++){
				fs>>temp;
				vdat[j]=temp;
				//cout<<temp<<endl;
			}
		//look at the type of particle and put it in the corresponding index, then do calculations e.g. pt
		
		vf->setTLVector((int)vdat[1],vdat[6],vdat[7],vdat[8],vdat[10]);
	
		}
	
	}
	//vf->printVectors();
	evtp=vf->arrPtr;

}
