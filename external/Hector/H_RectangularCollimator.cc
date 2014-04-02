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

/// \file H_RectangularCollimator.cc
/// \brief R-Collimators for the LHC beamline.

// local #includes
#include "H_RectangularCollimator.h"
#include "H_TransportMatrices.h"

void H_RectangularCollimator::init() {
	// must be in public section
		setTypeString();
		setMatrix(0,MP,QP);
	return;
}

H_RectangularCollimator::H_RectangularCollimator(const double s, const double l) :H_Drift(s,l){
	type = RCOLLIMATOR;
	init();
	return;
}

H_RectangularCollimator::H_RectangularCollimator(const string& nameE, const double s, const double l) :H_Drift(nameE,s,l){
	type = RCOLLIMATOR;
	init();
	return;
}

std::ostream& operator<< (std::ostream& os, const H_RectangularCollimator& el) {
	os << el.typestring << el.name << "\t\t at s = " << el.fs << "\t length = " << el.element_length <<endl;
	if(el.element_aperture->getType()!=NONE) {
        	os << *(el.element_aperture) << endl;
	}

	if(el.element_length<0 && VERBOSE) os <<"<H_RectangularCollimator> ERROR : Interpenetration of elements !"<<endl; 
	else if(el.element_length==0 && VERBOSE) os <<"<H_RectangularCollimator> WARNING : 0-length "<< el.name << " !" << endl;
   return os;
}

//void H_RectangularCollimator::setMatrix(const float eloss, const float p_mass, const float p_charge) {
void H_RectangularCollimator::setMatrix(const float , const float , const float ) {
	element_mat = driftmat(element_length);
	return ;
}

H_RectangularCollimator* H_RectangularCollimator::clone() const {
	H_RectangularCollimator* temp_coll = new H_RectangularCollimator(name,fs,element_length);
	temp_coll->setAperture(element_aperture);
	temp_coll->setX(xpos);
	temp_coll->setY(ypos);
	temp_coll->setTX(txpos);
	temp_coll->setTY(typos);
	temp_coll->setBetaX(betax);
	temp_coll->setBetaY(betay);
	return temp_coll;
}

