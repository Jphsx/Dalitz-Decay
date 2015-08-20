#include "vectorFactory.h"
using namespace std;
vectorFactory::vectorFactory(){}
vectorFactory::vectorFactory(double nu){
	gen = new eventGenerator(nu);
}
vectorFactory::vectorFactory(double nu, const char* CMfilepath){
	gen = new eventGenerator(nu,CMfilepath);
}
//reads a hepevt file from a specified path then builds the particle struct package and stores locally
void vectorFactory::readVector(const char* path){
	//this method can only be used to read a single event (it is redefined in minhelper) and is deprecated with N events
	
	ifstream f (path);
	
	int numPs;//number of particles in the event, as part of the hepevt formatting for each event
	double temp;
	//hardcoded sizes should be avoided
	double vdat[11];// array to store the pieces of the hepevt file on read
	
	

	if(f.is_open()){
	//used to set the number of iterations to read for each event
	f>>numPs;

	//instantiates a particle struct array with n+1 number of particles
	ParticleParameters *ptr = new ParticleParameters [numPs+1]; 
	//array struct instance is then assigned globally to the vector factory class
	arrPtr = ptr;
	
	//iterates through the events in from designated file and builds the struct
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
//helper function to build particle struct
void vectorFactory::setTLVector(int pID,double px,double py,double pz, double m){
	//conditionals look at the particle ID builds a single struct by setting the 4vector and passes the struct to make calculations
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
//performs some calculations to completely populate a particle struct i.e. the most commonly used variables
void vectorFactory::makeCalculations(ParticleParameters& particle){
	//sets transvers momentum, the polar displacement with respect to the z-axis(electron beam direction) in the simulated accelerator, and sets the measured curvature xm
	if(particle.pID==11 || particle.pID==-11){
		double Pt = TMath::Sqrt(particle.v.Px()*particle.v.Px() + particle.v.Py()*particle.v.Py());
		double p = TMath::Sqrt(Pt*Pt + particle.v.Pz()*particle.v.Pz());
		double thet = TMath::ACos(particle.v.Pz()/p);
	
		particle.pt=Pt;
		particle.theta=thet;
		particle.x_m=1/Pt;
	}
	//
	if(particle.pID==22){
		particle.pt=0.0;
		particle.theta=0.0;//this values needs to be checked and updated, however p3 theta is unused in all calculations
		particle.x_m=particle.v.E();
	}
	
}
// uses event generator to generate a particle event and builds a particle struct array
void vectorFactory::createEvent(double M, double m_e, double initial_p){
	//hardcoded to 3 particles
	ParticleParameters *ptr = new ParticleParameters [4];
	//have to manually set the PID for event generator (which builds the 4vectors instead of setTLVector)
	setPID(ptr);
	//event generator returns a 4vector array so one is instantiated locally to be returned to
	TLorentzVector* evtP = new TLorentzVector[4];
	//assigns the particle struct array globally
	arrPtr = ptr;
	//call generate event with the initial conditions, pion mass, electron mass, and initial momentum
	evtP = gen->generateEvent(M,m_e,initial_p);
	for(int i=1; i<=3; i++){
		//sets each 4 vector to its assigned position in the struct array
		arrPtr[i].v=evtP[i];
		//perfoms calculations each struct element 
		makeCalculations(arrPtr[i]);
	} 
}
//defines a particle ID values for a given particle struct array
void vectorFactory::setPID(ParticleParameters* evtP){
	for(int i=1; i<=3; i++){
		if(i==1) evtP[i].pID=-11;
		if(i==2) evtP[i].pID=11;
		if(i==3) evtP[i].pID=22;
	}
}
//creates a deep copy of struct a in struct b
void vectorFactory::copyObject(ParticleParameters* evtP_a, ParticleParameters* evtP_b){
	for(int i=1; i<=3; i++){
		evtP_b[i].v = evtP_a[i].v;
		evtP_b[i].theta = evtP_a[i].theta;
		evtP_b[i].pt = evtP_a[i].pt;
		evtP_b[i].x_m = evtP_a[i].x_m;
		evtP_b[i].pID = evtP_a[i].pID;
	}
}
/*int main(){
	double m_e=0.000511;
	double M=0.13497;
	
	vectorFactory *v = new vectorFactory(2*m_e/M);
	vectorFactory::ParticleParameters* evtP; 
	
	
	//v->readVector("dalitz_0mm_unsmeared.hepevt");
	v->createEvent(M,m_e,10.0);
	evtP= v->arrPtr;
	cout<<"about to print"<<endl;
	v->printVectors();
	TLorentzVector v0 = evtP[1].v+evtP[2].v+evtP[3].v;
	cout<<"V0 boosted " <<v0.Px()<< " "<<v0.Py()<< " "<<v0.Pz()<< " "<<v0.P()<< " "<<v0.E()<< " "<<v0.M()<< endl;
	
}*/
