/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_RectEllipticAperture.cc
/// \brief Defines the Rect-Elliptic aperture of beamline elements.

// C++ #includes
#include <iostream>

// C #includes
#include <cmath> // needed for fabs

// ROOT #includes
//#include "TEllipse.h"
//#include "TPave.h"

// local #includes
#include "H_RectEllipticAperture.h"
using namespace std;

H_RectEllipticAperture::H_RectEllipticAperture(const float x1, const float x2, const float x3, const float x4, const float posx, const float posy) :H_Aperture(RECTELLIPSE,((x1==0)?x3:x1),((x2==0)?x4:x2),x3,x4,posx,posy) {
	/// @param x1, x2, x3, x4 are the geometrical parameters of the rect-ellipse shape
	/// @param posx, posy defines the (x,y) of the center of the shape
	type= RECTELLIPSE;
}

void H_RectEllipticAperture::draw() const {
/*	TEllipse* te = new TEllipse(fx,fy,x3,x4);
	TPave* tp = new TPave(fx-x1,fy-x2,fx+x1,fy+x2,1);
        te->SetLineColor(2);
        te->SetFillColor(2);
	te->SetFillStyle(3004);
        te->Draw();
        tp->SetLineColor(2);
        tp->SetFillColor(2);
	tp->SetFillStyle(3005);
        tp->Draw();	
	return;
*/
}

bool H_RectEllipticAperture::isInside(const float x, const float y) const {
	return((fabs(fx-x)<x1)&&(fabs(fy-y)<x2)&&(((x-fx)/x3)*((x-fx)/x3)+((y-fy)/x4)*((y-fy)/x4)<1));
}

void H_RectEllipticAperture::printProperties() const {
	cout << "Aperture shape:" << getTypeString() << ", parameters " <<x1<<", "<<x2<<", "<<x3<<", "<<x4<< endl;
	cout << " \t Center : "<<fx<<", "<<fy<<endl;
	return;
}
