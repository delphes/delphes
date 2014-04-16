/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_RomanPot.cc
/// \brief Roman pot class

// local #includes
#include "H_RomanPot.h"
#include "H_RectangularAperture.h"
#include "H_TransportMatrices.h"

void H_RomanPot::init() {
	// must be in public section
		setTypeString();
		setMatrix(0,MP,QP);
	return;
}

H_RomanPot::H_RomanPot(const double s, const double app) :H_OpticalElement(RP,s,0.,RP_LENGTH){
		init();
		H_RectangularAperture* rapp = new H_RectangularAperture(app,RP_HEIGHT,0,0);
		setAperture(rapp);
}

H_RomanPot::H_RomanPot(const string nameE, const double s, const double app) :H_OpticalElement(nameE,RP,s,0.,RP_LENGTH){
		init();
		H_RectangularAperture* rapp = new H_RectangularAperture(app,RP_HEIGHT,0,0);
		setAperture(rapp);
}

void H_RomanPot::printProperties() const {
		cout << typestring << name;
		cout << "\t\t at s = " << fs;
                if(element_aperture->getType()!=NONE) {
                        cout <<"\t aperture type = " << element_aperture->getTypeString();
                        element_aperture->printProperties();
                }

		cout<<endl;
}

void H_RomanPot::setMatrix(const float eloss, const float p_mass, const float p_charge) const {
		*element_mat = driftmat(0);
	return ;
}
