#ifndef _MIN_HELPER
#define _MIN_HELPER
#include "../vectorFactory/vectorFactory.h"
//#include "../Chisq/DalitzChiSq.h"

//! Since MINUIT2 cannot be associated with a class when using an algorithm this class reads the .hepevt smeared and unsmeared files similar to vectorFactory. Then populates a particle struct array for use as global variables in the minimization implementation.
class MinHelper{
	public:
	/*! Particle struct array that is intended for the unsmeared event */
	vectorFactory::ParticleParameters* evtP_actual;
	/*! Particle struct array that is intended for the smeared event */
	vectorFactory::ParticleParameters* evtP;
	//f1 is smeared f2 is actual

	/*! Local global filestream f1 which is used to read the smeared event */
	ifstream f1;
	/*! Local global filestream f2 which is used to read the unsmeared event */
	ifstream f2;
	
	/*! Local instance of vectorFactory for use of struct population methods, each instance is associated with one filestream */
	vectorFactory* factory1;
	/*! Second instance of vectorFactory for use of struct population methods, each instance is associated with one filestream */
	vectorFactory* factory2;
	
	//! Method that populates a reference to the passed in struct with the associated filestream 
	/*!
		\param fs the filestream of the .hepevt file which is to be read from, the filestream to be passed in should be one of the classe's members
		\param vf the vectorFactory instance whose helper methods populate the given struct array
		\param evtp the reference to the address of the struct array that is to be populated
	*/
	void populateParticles(ifstream& fs,vectorFactory* vf,vectorFactory::ParticleParameters*& evtp);
	//! Constructor that sets the filestreams according to the filepaths passed in.  F1 should be the smeared event and F2 should be the unsmeared event
	MinHelper(const char* f1, const char* f2);
	
};
#endif
