#include "eventGenerator.h"
//constructor that initializes the random number generator and MonteCarlo xy generator with its necessary boundary conditions
eventGenerator::eventGenerator(double nu){
	RNG = new TRandom1();
	XYgenerator* temp = new XYgenerator(nu);
	//xygen = new XYgenerator(nu);
	xygen = temp;
	CMwrite=false;
}
eventGenerator::eventGenerator(double nu, const char* CMfilename){
	RNG = new TRandom1();
	XYgenerator* temp = new XYgenerator(nu);
	//xygen = new XYgenerator(nu);
	xygen = temp;
	CMwrite=true;
	CMfs.open(CMfilename);
}
// uses XY generator to get the x and y values then generates the three particles in the event, rotates the system about all three axes then boosts the particles into the lab frame
TLorentzVector* eventGenerator::generateEvent( double M, double m_e, double initial_p){
	double* xy = xygen->generateXYpair();
	// uses an array of four vectors to store the event temporarily
	TLorentzVector* evtP = new TLorentzVector[4];
	//iterates 3 times and makes each particle in the event individually
	for(int i=1; i<=3; i++){
		if(i == 1) evtP[i]=makePositron(xy[0],xy[1],M,m_e);
		if(i == 2) evtP[i]=makeElectron(xy[0],xy[1],M,m_e);
		if(i == 3) evtP[i]=makePhoton(xy[0],M);
	}
	/////////////////Must add function to write rotated CM events to output file////////////////////////////////
	if(CMwrite) writeEvent(CMfs,evtP);
	//rotate and boost
	rotateSystem(evtP);
	boostToLab(evtP, initial_p, M);
	//returns pointer to 4 vector array
	return evtP;
}
// creates a photon with a given x value and pion mass
TLorentzVector eventGenerator::makePhoton(double x, double M){
	//mass of gamma star is the virtual cumulative mass of the dielectric system before it is created
	double m_gamma_star = sqrt(x*M*M);
	//pcm is the scalar momentum of the photon in the CM frame which is dictated by the virtual particle mass and pion mass
	double pcm = (M*M - m_gamma_star*m_gamma_star)/(2*M);
	TLorentzVector v;
	//the default direction of the photon is the -x direction (for easy initial calculations, to be rotated later)
	v.SetXYZM(-pcm,0,0,0);
	//cout<<"PCM "<<pcm<<" "<<x<<endl;
	//cout<<v.E()<<" ENERGY"<<endl;
	return v;
}
// creates an electron with the given x and y values and appropriate masses
TLorentzVector eventGenerator::makeElectron(double x, double y,double M, double m_e){
	//calculates the 2 components of the electron 4 vector
	double P_lateral = (M/4) * ( (1-x) + y*(1+x) );
	double P_transverse = sqrt( ((x*M*M/4) * (1-y*y)) - m_e*m_e );
	//cout<<"Ps elec "<<P_lateral<<" "<< P_transverse <<" y "<<y<<" x "<<x<<endl;
	TLorentzVector v;
	//the components are set in the (x,y) plane which are not the same as the x,y passed in to calculate the components
	v.SetXYZM(P_lateral, P_transverse, 0, m_e);
	return v;
	
}
// creates a positron with the give x and y values and appropriate masses
TLorentzVector eventGenerator::makePositron(double x, double y, double M, double m_e){
	//calculates the 2 components of the electron 4 vector
	double P_lateral = (M/4) * ( (1-x) - y*(1+x) );
	//the transverse momentum of the positron is equal and opposite of electron transverse momentum
	double P_transverse = -sqrt( ((x*M*M/4) * (1-y*y)) - m_e*m_e);
	TLorentzVector v;
	//the components are set in the (x,y) plane which are not the same as the x,y passed in to calculate the components
	v.SetXYZM(P_lateral,P_transverse, 0, m_e);
	//cout<<"Ps posi "<<P_lateral<<" "<< P_transverse <<" y "<<y<<" x "<<x<<endl;
	return v; 
}
//rotates the 3 particle system about the x,y,z axes
void eventGenerator::rotateSystem(TLorentzVector*& evtP){
	//first rotates di-electron about x-axis
	double rads =2.0*M_PI*RNG->Uniform(0.0,1.0);
	for(int i=1; i<=2; i++){
		evtP[i].RotateX(rads);
	}
	//now rotate entire system about each axis
	rads =2.0*M_PI*RNG->Uniform(0.0,1.0);
	for(int i=1; i<=3; i++){
		evtP[i].RotateY(rads);
	}
	rads =2.0*M_PI*RNG->Uniform(0.0,1.0);
	for(int i=1; i<=3; i++){
		evtP[i].RotateZ(rads);
	}	
	//x-axis rotation is last because photon is originally fixed to the -x axis
	rads =2.0*M_PI*RNG->Uniform(0.0,1.0);
	for(int i=1; i<=3; i++){
		evtP[i].RotateX(rads);
	}
	//cout<<evtP[3].E()<<" rot energy"<<endl;
}
//boosts the CM particles into the lab frame with a given initial lab momentum
void eventGenerator::boostToLab(TLorentzVector*& evtP, double initial_p, double M){
	//magnitude of boost vector
	double beta = initial_p/sqrt(initial_p*initial_p + M*M);
	//random directions generated for the boost vector
	double phi = 2.0 * M_PI * RNG->Uniform(0.0,1.0);
	double costh = RNG->Uniform(-1.0,1.0);
	double theta = acos(costh);
	TVector3 b(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
	b.SetMag(beta);
	//iterate through event container and boost each particle into the lab
	for(int i=1; i<=3; i++){
		evtP[i].Boost(b);
	}
	//cout<<evtP[3].E()<<" boost engergy"<<endl;
	
}
//write the CM event to the CM .hepevt file
void eventGenerator::writeEvent(ostream& f1, TLorentzVector*& evtP){
		int quickPIDarr[4] ={-1,-11,11,22};
		f1<<3<<endl;
		for(int i=1; i<=3; i++){
			f1<< 1 <<" "<< quickPIDarr[i] << " " << 0 << " "<< 0 << " " << 0 << " "<< 0 << " ";
			//printVector(evtP[i].v);
			f1<<evtP[i].Px()<<" "<<evtP[i].Py()<<" "<<evtP[i].Pz()<<" "<<evtP[i].E()<<" "<<evtP[i].M()<<endl;
		}
} 
//debugging method that prints 4 vector contents
void eventGenerator::vectorPrinter(TLorentzVector* v){
	for(int i=1; i<=3; i++){
		cout<<v[i].Px()<<" "<<v[i].Py()<<" "<<v[i].Pz()<<" "<<v[i].P()<<" "<<v[i].E()<<" "<<v[i].M()<<endl;
	}
}
//int main(){
//cout<<"that"<<endl;}
