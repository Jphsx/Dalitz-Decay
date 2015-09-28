#include "DalitzSmear.h"
using namespace std;

DalitzSmear::DalitzSmear(){
	RNG=new TRandom1();
	 scaleParameterCP=1e-6;
	 scaleParameterNP=1e-4;
}
DalitzSmear::DalitzSmear(TLorentzVector v, int classifier){
	RNG=new TRandom1();
	pid = classifier;
	 scaleParameterCP=1e-6;
	 scaleParameterNP=1e-4;
	
	v_reg.SetXYZM(v.Px(),v.Py(),v.Pz(),v.M());
}
TLorentzVector DalitzSmear::SmearVector(){
	double ptsm = getPtsm(v_reg);
	v_sme.SetXYZM(ptsm*cos(mathUtility::getPhi(v_reg)), ptsm*sin(mathUtility::getPhi(v_reg)), ptsm*(1/tan(mathUtility::getTheta(v_reg))), v_reg.M());
	return v_sme;
}
TLorentzVector DalitzSmear::SmearVector(TLorentzVector v){
	double ptsm = getPtsm(v);
	
	v_reg.SetXYZM(v.Px(),v.Py(),v.Pz(),v.M());
	v_sme.SetXYZM(ptsm*cos(mathUtility::getPhi(v)), ptsm*sin(mathUtility::getPhi(v)), ptsm*(1/tan(mathUtility::getTheta(v))), v.M());
	
	return v_sme;
}
void DalitzSmear::setpid(int id){
	pid=id;
}
void DalitzSmear::setScaleParameterNP(double sigma){
	scaleParameterNP=sigma;
}
void DalitzSmear::setScaleParameterCP(double sigma){
	scaleParameterCP=sigma;
}
double DalitzSmear::getPtsm(TLorentzVector v){
	if(pid == 11 || pid == -11){
		return 1/ getGauss(v);	
	}
	if(pid == 22){
		return getGauss(v) * sin(mathUtility::getTheta(v));
	}
}
double DalitzSmear::getGauss(TLorentzVector v){
	if(pid == 11 || pid == -11){
		return RNG->Gaus(1/mathUtility::getPt(v), sqrt(mathUtility::getVariance(v,pid)));	
	}
	if(pid == 22){
		return RNG->Gaus(v.E(),sqrt(mathUtility::getVariance(v,pid)));
	}
}
TVector3 DalitzSmear::getUnitVector(double i, double j, double k){
	double mag = getVectorMagnitude(i,j,k);
	TVector3 u(i/mag, j/mag, k/mag);
	return u;
}
double DalitzSmear::getVectorMagnitude(double i, double j, double k){
	return sqrt(i*i + j*j + k*k);
}
//the directional smearing component of the system, this vector lies on a plane that is normal to the original generator level direction
TVector3 DalitzSmear::getOmega_T(TVector3 v_1, TVector3 v_n){
	double phi_t = RNG->Uniform(0,M_PI);
	TVector3 w_t = cos(phi_t)*v_1 + sin(phi_t)*v_n;
	return w_t;
}
//the vector passed in should be the unit vector that describes the direction of the generator level particle
//v1 is perpendicular to the original direction of the particle
//v1 is also in the same plane as w_t
TVector3 DalitzSmear::getV_1(TVector3 uVector){
	
	double thetaV1 = atan( 1/tan(uVector.Theta()) );
	//restrict v1 polar angle to [0,pi]
	if(thetaV1 < 0 ) thetaV1 = thetaV1 + M_PI;

	//if(pid==22) cout<<"U PHI "<<uVector.Phi()<<endl;
	double phiV1 = uVector.Phi() - mathUtility::safeAcos( -1 / (tan(uVector.Theta())*tan(thetaV1)) );

	//restrict v1 azimuth to [0,2pi]
	if( phiV1 < 0 ) phiV1 = phiV1 + 2*M_PI;
	
	
	
	TVector3 v1(sin(thetaV1)*cos(phiV1), sin(thetaV1)*sin(phiV1), cos(thetaV1));
	return v1;
}
// vector functions mostly for intermediate checks
double DalitzSmear::getScalarProduct(TVector3 v1, TVector3 v2){
	return v1.X()*v2.X() + v1.Y()*v2.Y() + v1.Z()*v2.Z();
}
TVector3 DalitzSmear::getVectorProduct(TVector3 v1, TVector3 v2){
	return v1.Cross(v2);
}
double DalitzSmear::getDpsi(TLorentzVector v){
	double r1 = RNG->Uniform(0,1);
	double negatives;
	//returns random variate of Rayleigh distribution with scale parameter sigma define in header
	if(pid==-11 || pid==11){
		if(scaleParameterCP*sqrt(-2.0*log(r1))<0.0){
			negatives = scaleParameterCP*sqrt(-2.0*log(r1));
			cout<<" (CP) R1 "<<r1<< "DPSI = "<<negatives<<endl;
		}
		return scaleParameterCP*sqrt( -2.0 * log(r1) );
	}
	if(pid==22){
		if(scaleParameterNP*sqrt(-2.0*log(r1))<0.0){
			negatives = scaleParameterNP*sqrt(-2.0*log(r1));
			cout<<" (NP) R1 "<<r1<< "DPSI = "<<negatives<<endl;
		}
		return scaleParameterNP*sqrt( -2.0 * log(r1) );
	}
}
TLorentzVector DalitzSmear::smearDirection(TLorentzVector v){
	//build from the bottom up
	TVector3 u = getUnitVector(v.Px(), v.Py(), v.Pz());
	
	//verify the magnitude is actually 1
	//if(getVectorMagnitude(u.X(),u.Y(),u.Z()) != 1.0)cout<<"u magnitude != 1, ||u|| = "<<getVectorMagnitude(u.X(),u.Y(),u.Z())<<endl;
	

	TVector3 v1 = getV_1(u);
	
	//verify that v1 is perpendicular to u
	//if(getScalarProduct(v1,u) != 0.0) cout<<"v1 is not perpendicular to u, scalar product: "<<getScalarProduct(v1,u)<<endl;
	
	TVector3 vn = getVectorProduct(u,v1);

	TVector3 wt = getOmega_T(v1,vn);
	
	//verify that wt is also perpendicular to u
	//if(getScalarProduct(wt,u) != 0.0) cout<<"wt is not perpendicular to u, scalar product: "<<getScalarProduct(wt,u)<<endl;
	
	double dpsi = getDpsi(v);

	//check if dpsi is larger than pi/2, if it is discard it and get a new value, 
	//in a few photon smearing cases the dot product between generator and detector four vectors was <0
	while(dpsi >= M_PI/2){
		dpsi = getDpsi(v);
	}
	
	
	TVector3 vsm = cos(dpsi)*u + sin(dpsi)*wt;
	//cout<<"cout works here"<<endl;
	while( cos( vsm.Dot(u)) < 0.0 ){
		cout<<" psi >= pi/2 " <<endl;
		dpsi = getDpsi(v);
		vsm = cos(dpsi)*u + sin(dpsi)*wt;
	}
	


	
	double mag = getVectorMagnitude(v.Px(),v.Py(),v.Pz());

	v_sme.SetXYZM(mag*vsm.X(),mag*vsm.Y(),mag*vsm.Z(),v.M());
	
	return v_sme;
}
///vector printing method for the test framework
/*void printVector(TLorentzVector v){
	cout<<v.Px()<<" "<<v.Py()<<" "<<v.Pz()<<" "<<v.P()<<" "<<v.E()<<" "<<v.M()<<endl;
}
int main(){
	//testing framework for magnitude// directional smearing

	//3 random events from generator level lab frame
	/*3
	1 -11 0 0 0 0 0.118038 1.83901 0.534643 1.91878 0.000511
	1 11 0 0 0 0 0.300449 4.51956 1.32268 4.7187 0.000511
	1 22 0 0 0 0 0.242968 3.24431 0.853262 3.36343 -4.21468e-08
	3
	1 -11 0 0 0 0 -0.951265 1.56249 -0.283607 1.85114 0.000511
	1 11 0 0 0 0 -0.224976 0.370414 -0.0654395 0.438296 0.000511
	1 22 0 0 0 0 -3.75269 6.61996 -1.24913 7.71148 1.19209e-07
	3
	1 -11 0 0 0 0 0.0838038 4.26598 -2.36552 4.87866 0.000511
	1 11 0 0 0 0 0.0734162 3.72949 -2.07093 4.26652 0.000511
	1 22 0 0 0 0 0.0119472 0.727345 -0.450668 0.855731 2.35608e-08
	*/
	//TLorentzVector test1_a,test1_b,test1_c, test2_a,test2_b,test2_c;
	//TLorentzVector test1_asm, test1_bsm, test1_csm, test2_asm, test2_bsm, test2_csm;
/*
	TLorentzVector* test1gen = new TLorentzVector[3]; //{particle 1,2,3 at index 0,1,2}
	TLorentzVector* test2gen = new TLorentzVector[3];
	TLorentzVector* test1sm = new TLorentzVector[3];
	TLorentzVector* test2sm = new TLorentzVector[3];

	test1gen[0].SetXYZM(0.118038, 1.83901, 0.534643, 0.000511);
	test1gen[1].SetXYZM(0.300449, 4.51956, 1.32268, 0.000511);
	test1gen[2].SetXYZM(0.242968, 3.24431, 0.853262, -4.21468e-08);

	test2gen[0].SetXYZM(-0.951265, 1.56249, -0.283607, 0.000511);
	test2gen[1].SetXYZM(-0.224976, 0.370414, -0.0654395, 0.000511);
	test2gen[2].SetXYZM(-3.75269, 6.61996, -1.24913, 1.19209e-07);

	DalitzSmear* sm = new DalitzSmear();
	//smear 1a
	sm->setpid(-11);
	test1sm[0] = sm->smearDirection(test1gen[0]);
	//smear 1b
	sm->setpid(11);
	test1sm[1] = sm->smearDirection(test1gen[1]);
	//smear 1c
	sm->setpid(22);
	test1sm[2] = sm->smearDirection(test1gen[2]);

	//now smear the second event
	sm->setpid(-11);
	test2sm[0] = sm->smearDirection(test2gen[0]);
	//smear 2b
	sm->setpid(11);
	test2sm[1] = sm->smearDirection(test2gen[1]);
	//smear 2c
	sm->setpid(22);
	test2sm[2] = sm->smearDirection(test2gen[2]);

	//display the results to screen
	cout<<"TEST #1"<<endl;
	cout<<setprecision(15);
	for(int i=0; i<3;i++){
		
		cout<<" particle: "<<i<<" generator level ";
		printVector(test1gen[i]);
		cout<<" THETA: "<<test1gen[i].Theta()<<" PHI: "<<test1gen[i].Phi()<<endl;
		cout<<endl;
		cout<<" particle: "<<i<<" detector level ";
		printVector(test1sm[i]);
		cout<<" THETA: "<<test1sm[i].Theta()<<" PHI: "<<test1sm[i].Phi()<<endl;
		cout<<endl;
		cout<<" theta - theta' "<<test1gen[i].Theta()- test1sm[i].Theta()<<endl;
		cout<<" phi - phi' " <<test1gen[i].Phi() - test1sm[i].Phi()<<endl;
		cout<<endl;

	}
	cout<<endl;
	cout<<"TEST #2"<<endl;
		for(int i=0; i<3;i++){
		
		cout<<" particle: "<<i<<" generator level ";
		printVector(test2gen[i]);
		cout<<" THETA: "<<test2gen[i].Theta()<<" PHI: "<<test2gen[i].Phi()<<endl;
		cout<<endl;
		cout<<" particle: "<<i<<" detector level ";
		printVector(test2sm[i]);
		cout<<" THETA: "<<test2sm[i].Theta()<<" PHI: "<<test2sm[i].Phi()<<endl;
		cout<<endl;
		cout<<" theta - theta' "<<test2gen[i].Theta()- test2sm[i].Theta()<<endl;
		cout<<" phi - phi' " <<test2gen[i].Phi() - test2sm[i].Phi()<<endl;
		cout<<endl;

	}
		

}*/
	
