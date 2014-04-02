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

/// \file H_Aperture.cc
/// \brief Classes aiming at simulating the aperture of the LHC beamline elements.

// C++ #includes
#include <iostream>

// local #includes
#include "H_Aperture.h"
using namespace std;

H_Aperture::H_Aperture() :
        type_(NONE), x1(0), x2(0), x3(0), x4(0), fx(0), fy(0),  aptypestring(getApertureString()) {
}

H_Aperture::H_Aperture(const int dtype, const float size1, const float size2, const float size3, const float size4, const float posx, const float posy) :
        type_(dtype), x1(size1), x2(size2), x3(size3), x4(size4), fx(posx), fy(posy), aptypestring(getApertureString()) {
    /// @param dtype defines the aperture shape
        /// @param size1, size2, size3, size4 are the geometrical sizes (length/width or great/small radii) in m
	/// @param posx, posy are the (x,y) coordinates of the center of the aperture [m]
}

H_Aperture::H_Aperture(const H_Aperture& ap) :
        type_(ap.type_), x1(ap.x1), x2(ap.x2), x3(ap.x3), x4(ap.x4), fx(ap.fx), fy(ap.fy), aptypestring(ap.aptypestring) {
}

H_Aperture& H_Aperture::operator=(const H_Aperture& ap) {
	if(this==&ap) return *this;
	type_ = ap.type_;
	x1 = ap.x1;
	x2 = ap.x2;
	x3 = ap.x3;
	x4 = ap.x4;
	fx = ap.fx;
	fy = ap.fy;
	aptypestring = ap.aptypestring;
	return *this;
}

std::ostream& operator<< (std::ostream& os, const H_Aperture& ap) {
	os << "Aperture shape:" << ap.aptypestring << ", parameters "<<ap.x1<<", "<<ap.x2<<", "<<ap.x3<<", "<<ap.x4<<endl;
	os << " \t Center : " << ap.fx << "," << ap.fy << endl;
	return os;
}

void H_Aperture::printProperties() const {
	cout << *this;
        return;
}

void H_Aperture::setPosition(const float posx, const float posy) {
	/// @param posx, posy are the (x,y) coordinates of the center of the apertures [m]
	// posx, posy, fx, fy = [m]
	fx = posx;
	fy = posy;
	return;
}

bool H_Aperture::isInside(const float , const float ) const {
        /// @param x, y are the (x,y) coordinates of the proton, in [m]
        //cout << "aperture::isInside" << endl;
        return false;
}


/*void H_Aperture::setApertureString() {
	switch (type_) {
		case NONE: aptypestring = NONENAME; break;
 		case RECTANGULAR: aptypestring = RECTANGULARNAME; break;
 		case ELLIPTIC: aptypestring = ELLIPTICNAME; break;
 		case CIRCULAR: aptypestring = CIRCULARNAME; break;
 		case RECTELLIPSE: aptypestring = RECTELLIPSENAME; break;
		default: aptypestring = NONENAME; break;
	}
}
*/

const string H_Aperture::getApertureString() const {
	string str;
        switch (type_) {
                case NONE: str = NONENAME; break;
                case RECTANGULAR: str = RECTANGULARNAME; break;
                case ELLIPTIC: str = ELLIPTICNAME; break;
                case CIRCULAR: str = CIRCULARNAME; break;
                case RECTELLIPSE: str = RECTELLIPSENAME; break;
                default: str = NONENAME; break;
        }
	return str;
}

