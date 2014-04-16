/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_EllipticAperture.cc
/// \brief Defines the elliptic aperture of beamline elements.

// C++ #includes
#include <iostream>

// ROOT #includes
//#include "TEllipse.h"

// local #includes
#include "H_Parameters.h"
#include "H_Aperture.h"
#include "H_EllipticAperture.h"
using namespace std;

H_EllipticAperture::H_EllipticAperture(const float x1, const float x2, const float posx, const float posy) :H_Aperture(ELLIPTIC,x1,x2,0,0,posx,posy) {
        /// @param x1, x2 are the length and height of the elliptic shape
        /// @param posx, posy are the (x,y) coordinates of the center of the ellipse
}

bool H_EllipticAperture::isInside(const float x, const float y) const {
        /// @param x, y are the (x,y) coordinates of the proton
	return (((x-fx)/x1)*((x-fx)/x1) + ((y-fy)/x2)*((y-fy)/x2) < 1);
}

void H_EllipticAperture::draw() const {
/* 	Ellipse* te = new TEllipse(fx,fy,x1,x2);
	te->SetFillStyle(3003);
	te->SetLineColor(2);
	te->SetFillColor(2);
	te->Draw();
	return;
*/
}

void H_EllipticAperture::printProperties() const {
	cout<< "Aperture shape:" << getTypeString() << ", ellipse axes : "<<x1<<", "<<x2<<endl;
	cout << " \t Center : " << fx << "," << fy << endl;
	return;
}
