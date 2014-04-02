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

H_RomanPot::H_RomanPot(const double s, const double app) :H_Drift(s,RP_LENGTH){
		type = RP;
		init();
		if(element_aperture) delete element_aperture;
		element_aperture = new H_RectangularAperture(app,RP_HEIGHT,0,0);
}

H_RomanPot::H_RomanPot(const string& nameE, const double s, const double app) :H_Drift(nameE,s,RP_LENGTH){
		type = RP;
		init();
		if(element_aperture) delete element_aperture;
		element_aperture = new H_RectangularAperture(app,RP_HEIGHT,0,0);
}

std::ostream& operator<< (std::ostream& os, const H_RomanPot& el) {
	os << el.typestring << el.name << "\t\t at s = " << el.fs;
	if(el.element_aperture->getType()!=NONE) {
		os << *(el.element_aperture) << endl;
	}
   return os;
}

void H_RomanPot::setMatrix(const float eloss, const float p_mass, const float p_charge) {
		element_mat = driftmat(0);
	return ;
}

H_RomanPot* H_RomanPot::clone() const {
	H_RomanPot* temp_rp = new H_RomanPot(name,fs,element_length);
	temp_rp->setAperture(element_aperture);
	temp_rp->setX(xpos);
	temp_rp->setY(ypos);
	temp_rp->setTX(txpos);
	temp_rp->setTY(typos);
	temp_rp->setBetaX(betax);
	temp_rp->setBetaY(betay);
	return temp_rp;
}

