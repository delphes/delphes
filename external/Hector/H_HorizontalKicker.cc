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

/// \file H_HorizontalKicker.cc
/// \brief Classes aiming at simulating horizontal kickers in beamline

// fk [rad] for kickers !!!!

// local #includes
#include "H_HorizontalKicker.h"
#include "H_TransportMatrices.h"

int kickers_on = 0;

void H_HorizontalKicker::setMatrix(const float eloss, const float p_mass, const float p_charge) {
	extern int kickers_on;
	if(kickers_on) {
		element_mat = hkickmat(element_length,fk,eloss,p_mass,p_charge);
	} else {
		element_mat = driftmat(element_length);
	}
	return ;
}

H_HorizontalKicker* H_HorizontalKicker::clone() const {
	H_HorizontalKicker* temp_kick = new H_HorizontalKicker(name,fs,fk,element_length);
	temp_kick->setAperture(element_aperture);
	temp_kick->setX(xpos);
	temp_kick->setY(ypos);
	temp_kick->setTX(txpos);
	temp_kick->setTY(typos);
	temp_kick->setBetaX(betax);
	temp_kick->setBetaY(betay);
	return temp_kick;
}

