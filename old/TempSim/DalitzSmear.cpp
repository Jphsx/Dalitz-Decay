#include "DalitzSmear.h"
using namespace std;

DalitzSmear::DalitzSmear(){RNG=new TRandom();}
DalitzSmear::DalitzSmear(TLorentzVector v, int classifier){
	RNG=new TRandom();
	pid = classifier;
	v_reg.SetXYZM(v.Px(),v.Py(),v.Pz(),v.M());
}
TLorentzVector DalitzSmear::SmearVector(){
	double ptsm = getPtsm(v_reg);
	v_sme.SetXYZM(ptsm*cos(getPhi(v_reg)), ptsm*sin(getPhi(v_reg)), ptsm*(1/tan(getTheta(v_reg))), v_reg.M());
	return v_sme;
}
TLorentzVector DalitzSmear::SmearVector(TLorentzVector v){
	double ptsm = getPtsm(v);
	
	v_reg.SetXYZM(v.Px(),v.Py(),v.Pz(),v.M());
	v_sme.SetXYZM(ptsm*cos(getPhi(v)), ptsm*sin(getPhi(v)), ptsm*(1/tan(getTheta(v))), v.M());
	
	return v_sme;
}
void DalitzSmear::setpid(int id){
	pid=id;
}
double DalitzSmear::getPt(TLorentzVector v){
	return sqrt(v.Px()*v.Px() + v.Py()*v.Py());
}
double DalitzSmear::getTheta(TLorentzVector v){
	return acos(v.Pz()/get_p(v));
}
double DalitzSmear::get_p(TLorentzVector v){
	return sqrt(getPt(v)*getPt(v) + v.Pz()*v.Pz());
}
double DalitzSmear::getPhi(TLorentzVector v){
	return atan2(v.Py(),v.Px());
}
double DalitzSmear::getPtsm(TLorentzVector v){
	if(pid == 11 || pid == -11){
		return 1/ getGauss(v);	
	}
	if(pid == 22){
		return getGauss(v) * sin(getTheta(v));
	}
}
double DalitzSmear::getGauss(TLorentzVector v){
	if(pid == 11 || pid == -11){
		return RNG->Gaus(1/getPt(v), sqrt(getVariance(v)));	
	}
	if(pid == 22){
		return RNG->Gaus(v.E(),sqrt(getVariance(v)));
	}
}
double DalitzSmear::getVariance(TLorentzVector v){
	if(pid == 11 || pid == -11){
		return (2e-5*2e-5) + pow((1e-3 * 1/getPt(v) * 1/sin(getTheta(v)) ),2);	
	}
	if(pid == 22){
		return pow(v.E(),2) * (pow(0.16/sqrt(v.E()),2) + (0.01*0.01));
	}
}
//int main(){}
	
