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

/// \file H_RectangularDipole.cc
/// \brief Classes aiming at simulating LHC beam rectangular dipoles

#include "H_RectangularDipole.h"
#include "H_TransportMatrices.h"

void H_RectangularDipole::setMatrix(const float eloss, const float p_mass, const float p_charge) {
	if (fk !=0 ) element_mat = rdipmat(element_length,fk,eloss,p_mass,p_charge);
	else  {
		element_mat = driftmat(element_length);
		if(VERBOSE) cout<<"<H_RectangularDipole> WARNING : k0= 0, drift-like dipole (" << name << ") !" << endl;
	}
	return ;
}

H_RectangularDipole* H_RectangularDipole::clone() const {
	H_RectangularDipole* temp_dip = new H_RectangularDipole(name,fs,fk,element_length);
	temp_dip->setAperture(element_aperture);
	temp_dip->setX(xpos);
	temp_dip->setY(ypos);
	temp_dip->setTX(txpos);
	temp_dip->setTY(typos);
	temp_dip->setBetaX(betax);
	temp_dip->setBetaY(betay);
	return temp_dip;
}

