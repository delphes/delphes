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

/// \file H_RectangularAperture.cc
/// \brief Defines the rectangular aperture of beamline elements.

// C++ #includes
#include <iostream>

// C #includes
#include <cmath> // needed for fabs

// ROOT #includes
#include "TPave.h"

// local #includes
#include "H_RectangularAperture.h"
using namespace std;

H_RectangularAperture::H_RectangularAperture(const float l, const float h, const float posx, const float posy) :H_Aperture(RECTANGULAR,l,h,0,0,posx,posy) {
	/// @param l, h are the length and height of the rectangular shape
	/// @param posx, posy are the (x,y) coordinates of the center of the rectangular shape
}

H_RectangularAperture* H_RectangularAperture::clone() const {
	return new H_RectangularAperture(x1,x2,fx,fy);
}

bool H_RectangularAperture::isInside(const float x, const float y) const {
	/// @param x, y are the (x,y) coordinates of the proton, in [m]
	return (fabs(x-fx)<x1&&fabs(y-fy)<x2);
}

void H_RectangularAperture::draw(const float scale) const {
	TPave* tp = new TPave((fx-x1)*scale,(fy-x2)*scale,(fx+x1)*scale,(fy+x2)*scale,1);
	tp->SetFillStyle(3003);
	tp->SetLineColor(2);
	tp->SetFillColor(2);
	tp->Draw();
	return;
}

std::ostream& operator<< (std::ostream& os, const H_RectangularAperture& ap) {
	os << "Aperture shape:" << ap.aptypestring << ", rectangle sides : " << ap. x1 <<", " << ap.x2 <<endl;
	os << " \t Center : " << ap.fx << "," << ap.fy << endl;
        return os;
}

