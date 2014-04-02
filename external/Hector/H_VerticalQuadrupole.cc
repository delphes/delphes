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

/// \file H_VerticalQuadrupole.cc
/// \brief Vertically focussing quadrupoles

// local #includes
#include "H_VerticalQuadrupole.h"

void H_VerticalQuadrupole::setMatrix(const float eloss, const float p_mass, const float p_charge) {
	if (fk<0)  { if(VERBOSE) cout<<"<H_VerticalQuadrupole> ERROR : k1 should be > 0 for H_VerticalQuadrupole (" << name << ")!"<<endl; }
		if (fk !=0 ) element_mat = vquadmat(element_length,fk,eloss, p_mass,p_charge);
		else  {
			element_mat = driftmat(element_length);
			if(VERBOSE) cout<<"<H_VerticalQuadrupole> WARNING : k1= 0, drift-like quadrupole (" << name << ") !" << endl;
   		}
	return ;
}

H_VerticalQuadrupole* H_VerticalQuadrupole::clone() const {
	H_VerticalQuadrupole* temp_quad = new H_VerticalQuadrupole(name,fs,fk,element_length);
	temp_quad->setAperture(element_aperture);
	temp_quad->setX(xpos);
	temp_quad->setY(ypos);
	temp_quad->setTX(txpos);
	temp_quad->setTY(typos);
	temp_quad->setBetaX(betax);
	temp_quad->setBetaY(betay);
	return temp_quad;
}
