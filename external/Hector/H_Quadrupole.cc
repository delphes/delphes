/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_Quadrupole.cc
/// \brief Classes aiming at simulating LHC beam quadrupoles.

#include "H_Quadrupole.h"

void H_Quadrupole::init() {
	// needed for call from H- and V-Quadrupoles constructor
	// must be in public section
		setTypeString();
		setMatrix(0,MP,QP);
	return;
}

void H_Quadrupole::printProperties() const{
	cout << typestring;
	cout << name;
	cout << "\t at s = " << fs;
	cout << "\t length = " << element_length;
	cout << "\t k1 = " << fk;
	cout<<endl;
        if(element_aperture->getType()!=NONE) {
        	cout <<"\t aperture type = " << element_aperture->getTypeString();
                element_aperture->printProperties();
	}

	if(element_length<0)  { if(VERBOSE) cout<<"\t ERROR : Interpenetration of elements !"<<endl; }
	if(element_length==0) { if(VERBOSE) cout<<"\t WARNING : 0-length "<< typestring << " !" << endl; }
}
