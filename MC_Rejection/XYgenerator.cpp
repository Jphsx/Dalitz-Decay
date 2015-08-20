#include "XYgenerator.h"
//constructor that sets the nu^2 bound globally and sets WMAX to a default value.  WMAX is based on the maximum value of the probability density function which was approximated outside of the program.  This also initializes the random number generator.
XYgenerator::XYgenerator(double bound){
	//temporary wmax
	WMAX = 9000;
	nu=bound;
	RNG = new TRandom1();
}
//returns the Beta value based on a given x randomly generated.  used to determine the boundary conditions for -Beta<y<Beta
double XYgenerator::getBeta(double x){
	return sqrt(1 - (nu*nu/x));
}
//retuns the value of the probability function given the parameters x and y
double XYgenerator::getProbability(double x, double y){
	double a= 0.033;
	return ( pow(abs(1+a*x),2)) *((pow(1-x,3))/(4*x)) *(1+y*y+ (nu*nu)/x );
}
//The MC rejection algorithm, returns true for a valid x,y set whose probability is greater than a random number generated between 0 and WMAX
bool XYgenerator::testProbability(double p, double x, double y){
	double test = RNG->Uniform(0.0,WMAX);
	if (getProbability(x,y) > WMAX){
		cout<<"EXCEEDED WMAX!!"<<endl;
	}
	if(getProbability(x,y) > test){
		return true;
	}
	else return false;
}
//creates an 2 element array of x,y values where nu^2<x<1 and -beta<y<beta, this function calls the MC rejection algorithm (test probability) and then returns the xy pair if they are valid according to the test algorithm
double* XYgenerator::generateXYpair(){
	static double xy[2];
	
	xy[0] = RNG->Uniform(nu*nu,1);
	xy[1] = RNG->Uniform(-getBeta(xy[0]),getBeta(xy[0]));
	
	//if(testProbability(getProbability(xy[0],xy[1]),xy[0],xy[1]))
	//	return xy;
	//else generateXYpair();
	while(!testProbability(getProbability(xy[0],xy[1]),xy[0],xy[1])){
		xy[0] = RNG->Uniform(nu*nu,1);
		xy[1] = RNG->Uniform(-getBeta(xy[0]),getBeta(xy[0]));
	}
	
	return xy;
}
/*int main(){ XYgenerator* test = new XYgenerator(0.007572053);
	double* xytest = test->generateXYpair();
	cout<<"VALUES: "<<xytest[0]<<" "<<xytest[1]<<endl;}
*/
