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

/// \file H_Kicker.cc
/// \brief Classes aiming at simulating kickers in LHC beamline.

// local #includes
#include "H_OpticalElement.h"
#include "H_Kicker.h"

void H_Kicker::init() {
	// needed for call from H- and V-Kickers constructor
	// must be in public section
		setTypeString();
		if (fk !=0 ) {
			setMatrix(0,MP,QP); 
		}
	return;
}

std::ostream& operator<< (std::ostream& os, const H_Kicker& el) {
	os << el.typestring << el.name <<"\t at s = "<< el.fs <<"\t length = "<< el.element_length;
	os<<"\t k0 = "<< el.fk <<endl;
	if(el.element_aperture->getType()!=NONE) {
		os << *(el.element_aperture) << endl;
        }

	if(el.element_length<0 && VERBOSE) os<<"<H_Kicker> ERROR : Interpenetration of elements !"<<endl; 
	else if(el.element_length==0 && VERBOSE) os<<"<H_Kicker> WARNING : 0-length "<< el.name << " !" << endl; 
  return os;
}
