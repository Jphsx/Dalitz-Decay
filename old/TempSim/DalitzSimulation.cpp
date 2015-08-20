#include "vectorFactory.h"
#include "DalitzSmear.h"
#include "DalitzChiSq.h"
#include "TLorentzVector.h"
#include <iostream>
#include <fstream>

void printVector(TLorentzVector v){
			cout<<v.Px()<<" "<<v.Py()<<" "<<v.Pz()<<" "<<v.E()<<" "<<v.M()<<endl;
	}
void updateObject(vectorFactory::ParticleParameters& particle){
	vectorFactory* v = new vectorFactory();
	v->makeCalculations(particle);
}
void printVectors(vectorFactory::ParticleParameters* arrPtr){
		for(int i=1; i<=3; i++){
			cout<<arrPtr[i].pID<<" "<<arrPtr[i].v.Px()<<" "<<arrPtr[i].v.Py()<<" "<<arrPtr[i].v.Pz()<<" "<<arrPtr[i].v.M()<<endl;
	cout<<arrPtr[i].pt<<" "<<arrPtr[i].theta<<" "<<arrPtr[i].x_m<<endl;
		}
	}
int main(){
	

	vectorFactory* generator = new vectorFactory();
	generator->readVector("../vectorFactory/dalitz_0mm_unsmeared.hepevt");
	vectorFactory* genCopy = new vectorFactory();
	genCopy->readVector("../vectorFactory/dalitz_0mm_unsmeared.hepevt");

	vectorFactory::ParticleParameters* evtP = generator->arrPtr;
	//dupe the orginal particle set for comparison later
	vectorFactory::ParticleParameters* evtP_actual = genCopy->arrPtr;

	DalitzSmear* corruptor = new DalitzSmear();
	for(int i=1;i<=3;i++){
		corruptor->setpid(evtP[i].pID);
		evtP[i].v = corruptor->SmearVector(evtP[i].v);
		updateObject(evtP[i]);
	}

	DalitzChiSq* chisq = new DalitzChiSq();
	chisq->generateContour(evtP_actual[1], evtP_actual[2], evtP[1], evtP[2], evtP[3],400, .001);
	
}
