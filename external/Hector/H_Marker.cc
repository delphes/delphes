  /* * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                         *
*                   --<--<--  A fast simulator --<--<--     *
*                 / --<--<--     of particle   --<--<--     *
*  ----HECTOR----<                                          *
*                 \ -->-->-- transport through -->-->--     *
*                   -->-->-- generic beamlines -->-->--     *
*                                                           *
* JINST 2:P09005 (2007)                                     *
*      X Rouby, J de Favereau, K Piotrzkowski (CP3)         *
*       http://www.fynu.ucl.ac.be/hector.html               *
*                                                           *
* Center for Cosmology, Particle Physics and Phenomenology  *
*              Universite catholique de Louvain             *
*                 Louvain-la-Neuve, Belgium                 *
 *                                                         *
   * * * * * * * * * * * * * * * * * * * * * * * * * * * */

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

std::ostream& operator<< (std::ostream& os, const H_Marker& el) {
	os << el.typestring << el.name << "\t\t at s = " << el.fs;
	if(el.element_aperture->getType()!=NONE) {
                os << *(el.element_aperture) << endl;
        }
   return os;
}

//void H_Marker::setMatrix(const float eloss, const float p_mass, const float p_charge) {
void H_Marker::setMatrix(const float , const float , const float ) {
		element_mat = driftmat(0);
	return ;
}

H_Marker* H_Marker::clone() const {
	H_Marker* temp_mkr = new H_Marker(name,fs);
	temp_mkr->setAperture(element_aperture);
	temp_mkr->setX(xpos);
	temp_mkr->setY(ypos);
	temp_mkr->setTX(txpos);
	temp_mkr->setTY(typos);
	temp_mkr->setBetaX(betax);
	temp_mkr->setBetaY(betay);
	return temp_mkr;
}
