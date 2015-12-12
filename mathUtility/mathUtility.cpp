#include "mathUtility.h"

//used by dalitzChiSq
double mathUtility::getVariance(pObject p){
	if(p.pID == 11 || p.pID == -11){
		//1/sin(theta) is due to the assumption imposed that the collider is cylindrical
		return (2e-5*2e-5) + pow((1e-3 * 1/p.pt * 1/sin(p.theta) ),2);	
	}
	if(p.pID == 22){
		return pow(p.x_m,2) * (pow(0.16/sqrt(p.x_m),2) + (0.01*0.01));
	}
}
//used by DalitzSmearing
double mathUtility::getVariance(TLorentzVector v, int pid){
	if(pid == 11 || pid == -11){
		return (2e-5*2e-5) + pow((1e-3 * 1/getPt(v) * 1/sin(getTheta(v)) ),2);	
	}
	if(pid == 22){
		return pow(v.E(),2) * (pow(0.16/sqrt(v.E()),2) + (0.01*0.01));
	}

}
//used by Numerical Minimization
double mathUtility::getVariance(int pID, double xm, double theta){
	if(pID == 11 || pID == -11){
		//1/sin(theta) is due to the assumption imposed that the collider is cylindrical
		return (2e-5*2e-5) + pow((1e-3 * xm * 1/sin(theta) ),2);	
	}
	if(pID == 22){
		return pow(xm,2) * (pow(0.16/sqrt(xm),2) + (0.01*0.01));
	}
}
double mathUtility::getZ12(pObject p1, pObject p2){
	double beta1 = p1.v.P()/p1.v.E();
	double beta2 = p2.v.P()/p2.v.E();
	return 2*((1/(beta1*beta2)) - getCosTheta(p1.v,p2.v));
}
//returns a function of opening angles between particle 2 and 3 used for the particle 3 constraint equation
double mathUtility::getZ23(pObject p2, pObject p3){
	double beta2 = p2.v.P()/p2.v.E();
	return 2*((1/beta2) - getCosTheta(p2.v,p3.v));
}
//returns a function of opening angles between particle 1 and 3 used for the particle 3 constraint equation
double mathUtility::getZ13(pObject p1, pObject p3){
	double beta1 = p1.v.P()/p1.v.E();
	return 2*((1/beta1) - getCosTheta(p1.v,p3.v));
}
double mathUtility::getTheta(TLorentzVector v){
	return acos(v.Pz()/get_p(v));
}
double mathUtility::getPhi(TLorentzVector v){
	return atan2(v.Py(),v.Px());
}
double mathUtility::getPt(TLorentzVector v){
	return sqrt(v.Px()*v.Px() + v.Py()*v.Py());
}
double mathUtility::get_p(TLorentzVector v){
	return sqrt(getPt(v)*getPt(v) + v.Pz()*v.Pz());
}
double mathUtility::getCosTheta(TLorentzVector v1, TLorentzVector v2){
	return (v1.Px()*v2.Px() + v1.Py()*v2.Py() + v1.Pz()*v2.Pz())/(v1.P()*v2.P());
}
double mathUtility::getX3constrained(double x1,double x2, pObject p1, pObject p2, pObject p3){
	double topterm =   M*M - 2*m_e*m_e - (getZ12(p1,p2)/(sin(p1.theta) * sin(p2.theta) *x1*x2));
	double bottomterm = ( (getZ13(p1,p3)/(sin(p1.theta) * x1)) + (getZ23(p2,p3)/(sin(p2.theta) * x2)) );
	return topterm/bottomterm;
}
double mathUtility::getX3constrainedMin2(pObject p1, pObject p2, pObject p3 ){
		TLorentzVector v12 = p1.v + p2.v;
		double Cosphi12 = getCosTheta(v12, p3.v);
		double topterm = M*M -v12.M()*v12.M();
		double bottomterm = 2*(v12.E() - v12.P()* Cosphi12);
	return topterm/bottomterm; 
		
}
//overloaded to calculate X3 with the supplied angle from the four vectors
double mathUtility::getX3constrainedMin2(TLorentzVector v12, TLorentzVector v3){
	double Cosphi12 = getCosTheta(v12,v3);
	double topterm = M*M - v12.M()*v12.M();
	double bottomterm = 2*(v12.E() - v12.P() * Cosphi12);
	return topterm/bottomterm;
}
//overloaded to calculate X3 with a passed in opening angle between v12 and v3
double mathUtility::getX3constrainedMin2(TLorentzVector v12, TLorentzVector v3, double psi12_3){
	
	double topterm = M*M - v12.M()*v12.M();
	double bottomterm = 2*(v12.E() - v12.P() * cos(psi12_3));
	return topterm/bottomterm;
}
//limits the domain to avoid floating point error
double mathUtility::safeAcos(double x){
	if (x < -1.0) x = -1.0 ;
	else if (x > 1.0) x = 1.0 ;
	return acos (x) ;
 }
double mathUtility::sign(double param){
	//std::cout<<"param :"<<param<<std::endl;	
	if(param<0.0) return -1;
	if(param>0.0) return 1;
	if(param == 0) return 0;
	
}
double mathUtility::bisection(double(*f)(double), double a, double b, double TOL, int N ){
	int n = 1;
	double sign_fa = sign( f(a) );

	double halfLength, p, sign_fp;
	
	while( n <= N ){
		 halfLength = (b-a)/2.0;
		 p = a + halfLength;
		//std::cout<<std::setprecision(15);
		//std::cout<<"p_"<<n<<" = "<< p<<"  hl: "<<halfLength<<std::endl;
		if(halfLength < TOL) return p;
		
		n++;
		sign_fp = sign( f(p) );
		if(sign_fp == 0.0) return p;
		if(sign_fp*sign_fa > 0) a=p;
		if(sign_fp*sign_fa < 0) b=p;
		
	}
	std::cout<<"iterated over all N N (may not have converged)="<<N<<std::endl;
	return p;
}
//int main(){}
