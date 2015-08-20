// Example on how to use the new Minimizer class in ROOT
//  Show usage with all the possible minimizers.
// Minimize the Rosenbrock function (a 2D -function) // Now X^2 of Dalitz Elimination
// This example is described also in
// http://root.cern.ch/drupal/content/numerical-minimization#multidim_minim
// input : minimizer name + algorithm name
// randomSeed: = <0 : fixed value: 0 random with seed 0; >0 random with given seed
//
//Author: L. Moneta Dec 2010
//Modified by Justin Anguiano 2015
#include "MinHelper.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <iostream>
#include "TMath.h"
#include <cmath>
#include "../vectorFactory/vectorFactory.cpp"
using namespace std;
//double z12=0.0;
//this one is explicit maybe it will work
double DalitzElimination(const double *vars){
		
	//positron : pt theta x1m : 2.59316 1.40359 0.38563
//electron : pt theta x1m : 2.81121 1.40353 0.355718
//z12 z23 z13 1.10663e-07 0.000730332 0.000741285
//photon energy 5.175
	
	double x3m=5.157;
	    double z12 = 1.10663e-07;
	    double z23 = 0.000730332;
	    double z13 = 0.000741285;
	double term1 = (pow(vars[0]-0.38563,2) ) / ( (2e-5*2e-5) + pow((1e-3 * 0.38563 * 1/TMath::Sin(1.40359) ),2) );
	double term2 = (pow(vars[1]-0.355718,2) )/ ( (2e-5*2e-5) + pow((1e-3 * 0.355718 * 1/TMath::Sin(1.40353) ),2) );	

	//now for e3
	double topterm =   0.13497*0.13497 - 2*0.000511*0.000511 - ( 4.07819e-08/(TMath::Sin(1.40359) * TMath::Sin(1.40353) *vars[0]*vars[1]));
	double bottomterm = ( (0.000741247/(TMath::Sin(1.40359) * vars[0])) + (0.0007303/(TMath::Sin(1.40353) * vars[1])) );
	double x3=topterm/bottomterm;
	
	double var3 = (x3m*x3m)*( (0.16/TMath::Sqrt(x3m))*(0.16/TMath::Sqrt(x3m)) + (0.01*0.01) );

	double term3 = pow(x3-x3m,2)/var3;
	

	return term1+term2+term3;

}

int NumericalMinimization(const char * minName = "Minuit2",
                          const char *algoName = "" ,
                          int randomSeed = -1)
{
   // create minimizer giving a name and a name (optionally) for the specific
   // algorithm
   // possible choices are:
   //     minName                  algoName
   // Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
   //  Minuit2                     Fumili2
   //  Fumili
   //  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS,
   //                              BFGS2, SteepestDescent
   //  GSLMultiFit
   //   GSLSimAn
   //   Genetic
   ROOT::Math::Minimizer* min =
      ROOT::Math::Factory::CreateMinimizer(minName, algoName);

   // set tolerance , etc...
   min->SetMaxFunctionCalls(100000000); // for Minuit/Minuit2
   min->SetMaxIterations(1000000);  // for GSL
   min->SetTolerance(0.0001);
   min->SetPrintLevel(1);

   // create funciton wrapper for minmizer
   // a IMultiGenFunction type
   ROOT::Math::Functor f(&DalitzElimination,2);
   double step[2] = {0.00001,0.00001};
   // starting point

   double variable[2] = { 0.380,0.350};
   //if (randomSeed >= 0) {
    //  TRandom2 r(randomSeed);
     // variable[0] = r.Uniform(-20,20);
      //variable[1] = r.Uniform(-20,20);
   //}

   min->SetFunction(f);

   // Set the free variables to be minimized!
   min->SetVariable(0,"x",variable[0], step[0]);
   min->SetVariable(1,"y",variable[1], step[1]);

   // do the minimization
   min->Minimize();

   const double *xs = min->X();
   std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): "
             << min->MinValue()  << std::endl;

   std::cout<<"cov mat :"<<min->CovMatrixStatus()<<" "<<min->CovMatrix(0,0)<<", "<<min->CovMatrix(0,1)<<", "<<min->CovMatrix(1,0)<<", "<<min->CovMatrix(1,1)<<endl;
   //cout::"x3 fit value	=	"
   // expected minimum is 0
   /*if ( min->MinValue()  < 1.E-4  && f(xs) < 1.E-4)
      std::cout << "Minimizer " << minName << " - " << algoName
                << "   converged to the right minimum" << std::endl;
   else {
      std::cout << "Minimizer " << minName << " - " << algoName
                << "   failed to converge !!!" << std::endl;
      Error("NumericalMinimization","fail to converge");
   }*/

   return 0;
}
int main(){
//positron : pt theta x1m : 2.59316 1.40359 0.38563
//electron : pt theta x1m : 2.81121 1.40353 0.355718
//z12 z23 z13 1.10663e-07 0.000730332 0.000741285
	//z12=1.10663e-07;
	DalitzChiSq* calcHelp = new DalitzChiSq();
	MinHelper* h = new MinHelper("../EventOutputs/DalitzSmearVectors.hepevt","../EventOutputs/DalitzActualVectors.hepevt");
	cout<<"!@#!@#!@#!@#!#!@#!@#!@#test"<<endl;
	cout<<h->evtP[1].v.Px()<<endl;
	cout<<calcHelp->getZ12(h->evtP[1],h->evtP[2])<<endl;
	
	const char * minName = "Minuit2";
	const char *algoName = "DalitzElimination";
	int randomSeed = -1;
	NumericalMinimization( minName ,
                          algoName  ,
                          randomSeed );
}

 
