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

/// \file H_Drift.cc
/// \brief Class aiming at simulating LHC beam drift.

// local includes
#include "H_Drift.h"
#include "H_TransportMatrices.h"


void H_Drift::init() {
	// must be in public section
	element_mat.ResizeTo(MDIM,MDIM);
	setTypeString();
	setMatrix(0,MP,QP);
	return;
}

std::ostream& operator<< (std::ostream& os, const H_Drift& el) {
	os << el.typestring << el.name << "\t\t at s = " << el.fs << "\t length = " << el.element_length <<endl;
	if(el.element_aperture->getType()!=NONE) {
		os << *(el.element_aperture) << endl;
        }
	if(el.element_length<0 && VERBOSE) os <<"<H_Drift> ERROR : Interpenetration of elements !"<<endl;
	else if(el.element_length==0 && VERBOSE) os <<"<H_Drift> WARNING : 0-length "<< el.name << " !" << endl;
   return os;
}

//void H_Drift::setMatrix(const float eloss, const float p_mass, const float p_charge) {
void H_Drift::setMatrix(const float , const float , const float ) {
	element_mat = driftmat(element_length);
	return ;
}

H_Drift* H_Drift::clone() const {
	H_Drift* temp_drift = new H_Drift(name,fs,element_length);
	temp_drift->setX(xpos);
	temp_drift->setY(ypos);
	temp_drift->setTX(txpos);
	temp_drift->setTY(typos);
	temp_drift->setBetaX(betax);
	temp_drift->setBetaY(betay); 
	return temp_drift;
}

