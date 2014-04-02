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

/// \file H_Dipole.cc
/// \brief Class aiming at simulating LHC beam dipoles

//local includes
#include "H_Dipole.h"

void H_Dipole::init() {
	// needed for call from R- and S-Dipoles constructor
	// must be in public section
	    element_mat.ResizeTo(MDIM,MDIM);
		setTypeString();
		if (fk !=0 ) setMatrix(0,MP,QP);
	return;
}

std::ostream& operator<< (std::ostream& os, const H_Dipole& el) {
	os << el.typestring << el.name <<"\t at s = "<< el.fs <<"\t length = "<< el.element_length <<"\t k0 = "<<el.fk << endl;
        if(el.element_aperture->getType()!=NONE) {
                os << *(el.element_aperture) << endl;
	}

	if(el.element_length<0)  { if(VERBOSE) os <<"<H_Dipole> ERROR : Interpenetration of elements !"<<endl; }
	if(el.element_length==0) { if(VERBOSE) os <<"<H_Dipole> WARNING : 0-length "<< el.name << " !" << endl; }
   return os;
}
