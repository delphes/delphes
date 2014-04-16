/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_Marker.cc
/// \brief Class defining a marker in the beamline (e.g. for interaction point)

// local includes
#include "H_Marker.h"
#include "H_TransportMatrices.h"

void H_Marker::init() {
	// must be in public section
	setTypeString();
	setMatrix(0,MP,QP);
	return;
}

void H_Marker::printProperties() const {
	cout << typestring << name;
	cout << "\t\t at s = " << fs << endl;
	if(element_aperture->getType()!=NONE) {
		cout <<"\t aperture type = " << element_aperture->getTypeString();
                element_aperture->printProperties();
        }
}

void H_Marker::setMatrix(const float eloss, const float p_mass, const float p_charge) const {
		*element_mat = driftmat(0);
	return ;
}
