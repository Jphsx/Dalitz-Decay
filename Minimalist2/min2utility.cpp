
#include "min2utility.h"

//sigma_psi12_3 must be set for every time program is run
min2utility::min2utility(double M,TLorentzVector v12, TLorentzVector v3){
	//scaleParameterCP = sigma_sp;
	//this value was calculated from the simulation in scriptTools/Sigma_Psi123/Psi123RMS.cpp
	//sigma_psi12_3 = 9.96779238e-05;
	Mpi=M;
	x12_u = v12;
	x3_u = v3;
}
//secondary constructor, this needs you to use setvectors
min2utility::min2utility(double M){
	
	//sigma_psi12_3 = 9.96779238e-05;
	Mpi=M;
}
void min2utility::setSigma_psi12_3(double sigma){
	sigma_psi12_3 = sigma;
}
void min2utility::setVectors(TLorentzVector v12, TLorentzVector v3){
	x12_u = v12;
	x3_u = v3;
}
double min2utility::getChiSqMin2(double Mpi, TLorentzVector v12, TLorentzVector v3, double psi123, double psi123m){
	
	 //edit constrained to take optimal psivalue
	//std::cout<<"e3c "<<mathUtility::getX3constrainedMin2(v12,v3,psi123)<<endl;
	//std::cout<<"e3stdev "<<sqrt(mathUtility::getVariance(v3,22))<<endl;
	//std::cout<<"e3var "<<mathUtility::getVariance(v3,22)<<endl;
	
	return pow( (mathUtility::getX3constrainedMin2(v12,v3,psi123) - v3.E())/(sqrt(mathUtility::getVariance(v3,22))),2) + pow( (psi123 - psi123m)/(sigma_psi12_3),2 );
}
//the deriviative of the X^2 equation with respect to psi12,3 (the value to be adjusted) this function is to be passed into bisection to find the 
//optimal psi for X^2 minimization
double min2utility::Min2ChiPrime(double psi123){
	
	
	double E3cMin2 = mathUtility::getX3constrainedMin2(x12_u, x3_u, psi123);
	double firstTerm = -x12_u.P()*(Mpi*Mpi - x12_u.M()*x12_u.M())*sin(psi123);
	double secondTerm = (E3cMin2 - x3_u.E());
	double thirdTerm = mathUtility::getVariance(x3_u,22) * pow( (x12_u.E()- x12_u.P()*cos(psi123)), 2);
	
	double fourthTerm = (2*(psi123 - mathUtility::safeAcos(mathUtility::getCosTheta(x12_u,x3_u)))) / pow(sigma_psi12_3,2);
	return (firstTerm*secondTerm)/thirdTerm + fourthTerm;
}
double min2utility::sign(double param){
	//std::cout<<"param :"<<param<<std::endl;	
	if(param<0.0) return -1;
	if(param>0.0) return 1;
	if(param == 0) return 0;
	
}
//double min2utility::bisection(double min2utility::*f(double), double a, double b, double TOL, int N ){
//has some extra parameters for debugging
//double min2utility::bisection(double a, double b, double TOL, int N,
//TLorentzVector v12, TLorentzVector v3, double psi123m
// ){
double min2utility::bisection(double a, double b, double TOL, int N){	
	int n = 1;
	double sign_fa = sign( Min2ChiPrime(a) );

	double halfLength, p, sign_fp;
	
	while( n <= N ){ 	
		cout<<setprecision(9);
		//cout<< "psi123_"<<n<<" "<<p<<endl;
		 halfLength = (b-a)/2.0;
		 p = a + halfLength;
		//std::cout<<std::setprecision(15);
		//std::cout<<"p_"<<n<<" = "<< p<<"  hl: "<<halfLength<<std::endl;
		//std::cout<< getChiSqMin2(0.13497, v12, v3, p, psi123m)<<endl<<endl;

		if(halfLength < TOL) return p;
		
		n++;
		sign_fp = sign( Min2ChiPrime(p) );
		if(sign_fp == 0.0) return p;
		if(sign_fp*sign_fa > 0) a=p;
		if(sign_fp*sign_fa < 0) b=p;
		
	}
	std::cout<<"iterated over all N N (may not have converged)="<<N<<std::endl;
	return p;
}
//minimizes with the values supplied from the constructor for this particular instance of the class
//returns the optimal psi12,3 such that X^2 is minimized
double min2utility::MinimizeMin2(){
	//we can search between 0 and pi because the opening angle is boundable on this interval
	
	return bisection(0,M_PI,1e-15,10000);

}
/*int main(){
//the helper test framework, will look at 1 10GEV event, (10 GEV event #2) Actual & Smeared and test X^2 values and 
	TLorentzVector v12, v1,v2,v3, v12sm, v1sm, v2sm, v3sm;
	TFile *f = new TFile("min2test.root","RECREATE");

	v1.SetXYZM(-0.951264672, 1.56248772, -0.283606968, 0.000511);
	v2.SetXYZM(-0.224975998, 0.370413795, -0.0654394757, 0.000511);
	v3.SetXYZM( -3.75269207, 6.61996006, -1.24913409, 1.1920929e-07);
	
	//1 -11 0 0 0 0 -0.950990729 1.56203603 -0.283521362 1.85060114 0.000510999999
	//1 11 0 0 0 0 -0.224939823 0.370354192 -0.0654292572 0.438225285 0.000511
	//1 22 0 0 0 0 -3.99826618 7.04448659 -1.32873886 8.20831717 1.68587394e-07
	
	v1sm.SetXYZM(-0.950990729, 1.56203603, -0.283521362 , 0.000510999999);
	v2sm.SetXYZM(-0.224939823, 0.370354192, -0.0654292572 , 0.000511);
	v3sm.SetXYZM(-3.99826618, 7.04448659, -1.32873886 , 1.68587394e-07);

	v12 = v1+v2;
	v12sm = v1sm + v2sm;

	min2utility* test = new min2utility(0.13497, v12sm, v3sm);
    
        TH1D *hmin = new TH1D("hmin","hmin",31400,0.0,3.14);
   
	//plot the derivative function
	//TCanvas *c1 = new TCanvas("c1","x2prime",200,10,800,600);
	int n= 31400;
	double x[n],y[n];
	
	for(int i=0; i<n; i++){
		x[i] = double(i)*0.0001 + 0.00005;
		y[i] = test->Min2ChiPrime(x[i]);
		//cout<< x[i]<< " "<<y[i]<<endl;
                hmin->Fill(x[i],y[i]);		
	}
	//TMultiGraph *mg = new TMultiGraph();
	//TGraph* hchiprime = new TGraph(n,x,y);
	//hchiprime->Draw("AC*");
	//mg->Add(hchiprime);
	//mg->Draw("ACP");
        hmin->Draw();
	//c1->Print("X2prime.pdf");

        f->Write();
  	
	//double bestguessPSI = test->MinimizeMin2();
	double realPSI = mathUtility::safeAcos(mathUtility::getCosTheta(v12,  v3));
	double measuredPSI = mathUtility::safeAcos(mathUtility::getCosTheta(v12sm,  v3sm));

	double psidebug = test->bisection(0,M_PI,1e-15,10000, v12sm, v3sm, measuredPSI );

	cout<<setprecision(9);
	//cout<<"bisection Psi: "<<bestguessPSI<<" real Psi: "<<realPSI<<endl;
	//cout<<"X^2 minimization: "<< test->getChiSqMin2(0.13497, v12sm, v3sm, bestguessPSI, measuredPSI);
	//cout<<" X^2 real psi: "<< test->getChiSqMin2(0.13497, v12sm, v3sm, realPSI, measuredPSI)<<endl;
	
	//f->Write();

}*/

