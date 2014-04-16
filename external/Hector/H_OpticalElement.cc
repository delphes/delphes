/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_OpticalElement.cc
/// \brief Class aiming at describing any beam optical element.

// c++ #includes
#include <iostream>
#include <string>

// ROOT #includes
//#include "TPaveLabel.h"

// local #includes
#include "H_TransportMatrices.h"
#include "H_Aperture.h"
#include "H_OpticalElement.h"
using namespace std;

/// called by the constructors
void H_OpticalElement::init(const string nameE, const int typeE, const double s, const double k, const double l, H_Aperture* the_app) {
	// this is called by the constructors
	// must be in public section !
	// xpos and ypos are vectors with n point. They define the aperture shape of the optical element.
	name = nameE;
	fs = s;
 	fk = k;
	xpos = 0;
	ypos = 0;
	txpos = 0;
	typos = 0;
 	element_length = l;
 	type = typeE;
 	element_mat = new TMatrix(MDIM,MDIM);
	setAperture(the_app);
	
	if(element_length<0)  { if(VERBOSE) cout<<"\t ERROR : Interpenetration of elements !"<<endl; }
	if(element_length==0) { if(VERBOSE) cout<<"\t WARNING : 0-length element ! (" << name << ") " << " at " << fs << endl; }
	betax =0;
	betay =0;
}

H_OpticalElement::H_OpticalElement() {
	H_Aperture* the_app = new H_Aperture();
	init("",DRIFT,0.,0.,0.1,the_app);
}

H_OpticalElement::H_OpticalElement(const string nameE, const int typeE, const double s, const double k, const double l) {
	H_Aperture* the_app = new H_Aperture();
	init(nameE,typeE,s,k,l,the_app);
}

H_OpticalElement::H_OpticalElement(const string nameE, const int typeE, const double s, const double k, const double l, H_Aperture* the_app) {
	init(nameE,typeE,s,k,l,the_app);
}

H_OpticalElement::H_OpticalElement(const int typeE, const double s, const double k, const double l, H_Aperture* the_app) {
	init("",typeE,s,k,l,the_app);
}

H_OpticalElement::H_OpticalElement(const int typeE, const double s, const double k, const double l) {
	H_Aperture* the_app = new H_Aperture();
	init("",typeE,s,k,l,the_app);
}

H_OpticalElement::H_OpticalElement(const H_OpticalElement& el) {
	fs = el.fs;
	element_length = el.element_length;
	fk = el.fk;
	xpos = el.xpos;
	ypos = el.ypos;
	txpos = el.txpos;
	typos = el.typos;
	betax = el.betax;
	betay = el.betay;
	type = el.type;
	name = el.name;
	typestring = el.typestring;
	element_mat = new TMatrix(*(el.element_mat));
	element_aperture = new H_Aperture(*(el.element_aperture));
}

H_OpticalElement& H_OpticalElement::operator=(const H_OpticalElement& el) {
	if(this==&el) return *this;
        fs = el.fs;
        element_length = el.element_length;
        fk = el.fk;
        xpos = el.xpos;
        ypos = el.ypos;
		txpos = el.txpos;
		typos = el.typos;
        betax = el.betax;
        betay = el.betay;
        type = el.type;
        name = el.name;
        typestring = el.typestring;
	delete element_mat;
	delete element_aperture;
        element_mat = new TMatrix(*(el.element_mat));
        element_aperture = new H_Aperture(*(el.element_aperture));
	return *this;
}

void H_OpticalElement::setAperture(H_Aperture * ap) {
//	element_aperture = const_cast<H_Aperture*>(ap);
	element_aperture = ap;
	return;
}

void H_OpticalElement::printProperties() const {
		cout << typestring;
		cout << name;
		cout <<"\t at s = " << fs;
		cout <<"\t length = "<< element_length;
		if(fk!=0) cout <<"\t strength = " << fk;
		if(element_aperture->getType()!=NONE) { 
			cout <<"\t aperture type = " << element_aperture->getTypeString();
			element_aperture->printProperties();			
		}
	
		cout<<endl;
		if(element_length<0)  { if(VERBOSE) cout<<"\t ERROR : Interpenetration of elements !"<<endl; }
		if(element_length==0) { if(VERBOSE) cout<<"\t WARNING : 0-length "<< typestring << " !" << endl; }

 return;
}

void H_OpticalElement::showMatrix() const {
		printMatrix(element_mat);
	return;
}

void H_OpticalElement::drawAperture() const {
	element_aperture->draw();
	return;
}

TMatrix H_OpticalElement::getMatrix() const {
	return *element_mat;
}

TMatrix H_OpticalElement::getMatrix(const float eloss, const float p_mass, const float p_charge) const {
	setMatrix(eloss,p_mass,p_charge);
	return *element_mat;
}

void H_OpticalElement::draw(const float meight, const float height) const{
/*	/// @param meight is the minimal extend of the graph
	/// @param height is the maximal extend of the graph
	float x1 = getS();
	float x2 = getS() + getLength();
	float y1 = meight;
	float y2 = height;
	TPaveLabel* cur_box = new TPaveLabel(x1,y1,x2,y2,"");
	cur_box->SetBorderSize(1);
	cur_box->SetFillStyle(1001);
	cur_box->SetFillColor((int)getType());
	cur_box->Draw();
*/
}

