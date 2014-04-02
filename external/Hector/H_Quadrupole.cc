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

/// \file H_Quadrupole.cc
/// \brief Classes aiming at simulating LHC beam quadrupoles.

#include "H_Quadrupole.h"

void H_Quadrupole::init() {
	// needed for call from H- and V-Quadrupoles constructor
	// must be in public section
		element_mat.ResizeTo(MDIM,MDIM);
		setTypeString();
		setMatrix(0,MP,QP);
	return;
}

std::ostream& operator<< (std::ostream& os, const H_Quadrupole& el) {
	os << el.typestring  << el.name  << "\t at s = " << el.fs << "\t length = " << el.element_length << "\t k1 = " << el.fk <<endl;
        if(el.element_aperture->getType()!=NONE) {
        	os << *(el.element_aperture) << endl;
	}

	if(el.element_length<0 && VERBOSE) os<<"<H_Quadrupole> ERROR : Interpenetration of elements !"<<endl; 
	else if(el.element_length==0 && VERBOSE) os<<"<H_Quadrupole> WARNING : 0-length "<< el.typestring << " !" << endl; 
   return os;
}
