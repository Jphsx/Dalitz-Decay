#include "DalitzSmear.h"
using namespace std;

DalitzSmear::DalitzSmear(){RNG=new TRandom1();}
DalitzSmear::DalitzSmear(TLorentzVector v, int classifier){
	RNG=new TRandom1();
	pid = classifier;
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
//int main(){}
	
