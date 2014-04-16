/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_RectangularCollimator.cc
/// \brief R-Collimators for the LHC beamline.

// local #includes
#include "H_RectangularCollimator.h"
#include "H_TransportMatrices.h"

void H_RectangularCollimator::init() {
	// must be in public section
		setTypeString();
		setMatrix(0,MP,QP);
	return;
}

H_RectangularCollimator::H_RectangularCollimator(const double s, const double l) :H_OpticalElement(RCOLLIMATOR,s,0.,l){
		init();
	return;
}

H_RectangularCollimator::H_RectangularCollimator(const string nameE, const double s, const double l) :H_OpticalElement(nameE,RCOLLIMATOR,s,0.,l){
		init();
	return;
}

void H_RectangularCollimator::printProperties() const {
	cout << typestring << name;
	cout << "\t\t at s = " << fs;
	cout << "\t length = " << element_length;
	cout<<endl;
	if(element_aperture->getType()!=NONE) {
        	cout <<"\t aperture type = " << element_aperture->getTypeString();
                	element_aperture->printProperties();
	}

	if(element_length<0)  { if(VERBOSE) cout<<"\t ERROR : Interpenetration of elements !"<<endl; }
	if(element_length==0) { if(VERBOSE) cout<<"\t WARNING : 0-length "<< name << " !" << endl; }
}

void H_RectangularCollimator::setMatrix(const float eloss, const float p_mass, const float p_charge) const {
	//	cout<<"\t WARNING : H_RectangularCollimator matrices not yet implemented ! " << endl;
		*element_mat = driftmat(element_length);
	return ;
}
