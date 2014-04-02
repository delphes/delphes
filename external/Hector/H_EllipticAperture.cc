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

/// \file H_EllipticAperture.cc
/// \brief Defines the elliptic aperture of beamline elements.

// C++ #includes
#include <iostream>

// ROOT #includes
#include "TEllipse.h"

// local #includes
#include "H_Parameters.h"
#include "H_Aperture.h"
#include "H_EllipticAperture.h"
using namespace std;

H_EllipticAperture::H_EllipticAperture(const int type, const float l, const float h, const float posx, const float posy) :
	H_Aperture(type,l,h,0,0,posx,posy) {
        /// @param type is the aperture type (ELLIPTIC or CIRCULAR)
        /// @param l, h are the length and height of the elliptic shape
        /// @param posx, posy are the (x,y) coordinates of the center of the ellipse

	if (type!= ELLIPTIC && type != CIRCULAR) { 
		cout << "Warning: trying to define an EllipticalAperture which is neither elliptical nor circular." << endl;
		cout << "'Elliptical' type forced\n"; 
		type_ = ELLIPTIC; 
		aptypestring = getApertureString(); 
	}
}

H_EllipticAperture::H_EllipticAperture(const float l, const float h, const float posx, const float posy) :
	H_Aperture(ELLIPTIC,l,h,0,0,posx,posy) {
        /// @param l, h are the length and height of the elliptic shape
        /// @param posx, posy are the (x,y) coordinates of the center of the ellipse
}

H_EllipticAperture* H_EllipticAperture::clone() const {
	return new H_EllipticAperture(x1,x2,fx,fy);
}

bool H_EllipticAperture::isInside(const float x, const float y) const {
        /// @param x, y are the (x,y) coordinates of the proton
	return (((x-fx)/x1)*((x-fx)/x1) + ((y-fy)/x2)*((y-fy)/x2) < 1);
}

void H_EllipticAperture::draw(const float scale) const {
	TEllipse* te = new TEllipse(fx*scale,fy*scale,x1*scale,x2*scale);
	te->SetFillStyle(3003);
	te->SetLineColor(39);
	te->SetFillColor(39);
	te->Draw("f");
	return;
}

std::ostream& operator<< (std::ostream& os, const H_EllipticAperture& ap) {
	os<< "Aperture shape:" << ap.aptypestring << ", ellipse axes : "<< ap.x1 <<", " << ap.x2 << endl;
	os << " \t Center : "  << ap.fx << "," << ap.fy << endl;
	return os;
}
