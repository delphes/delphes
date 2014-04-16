/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_RectangularAperture.cc
/// \brief Defines the rectangular aperture of beamline elements.

// C++ #includes
#include <iostream>

// C #includes
#include <cmath> // needed for fabs

// ROOT #includes
//#include "TPave.h"

// local #includes
#include "H_RectangularAperture.h"
using namespace std;

H_RectangularAperture::H_RectangularAperture(const float x1, const float x2, const float posx, const float posy) :H_Aperture(RECTANGULAR,x1,x2,0,0,posx,posy) {
	/// @param x1, x2 are the length and height of the rectangular shape
	/// @param posx, posy are the (x,y) coordinates of the center of the rectangular shape
}

bool H_RectangularAperture::isInside(const float x, const float y) const {
	/// @param x, y are the (x,y) coordinates of the proton, in [m]
	return (fabs(x-fx)<x1&&fabs(y-fy)<x2);
}

void H_RectangularAperture::draw() const {
/*	TPave* tp = new TPave(fx-x1,fy-x2,fx+x1,fy+x2,1);
	tp->SetFillStyle(3003);
	tp->SetLineColor(2);
	tp->SetFillColor(2);
	tp->Draw();
	return;
*/
}
void H_RectangularAperture::printProperties() const {
	cout << "Aperture shape:" << getTypeString() << ", rectangle Sides : "<<x1<<", "<<x2<<endl;
	cout << " \t Center : " << fx << "," << fy << endl;
	return;
}
