/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_Aperture.cc
/// \brief Classes aiming at simulating the aperture of the LHC beamline elements.

// C++ #includes
#include <iostream>

// local #includes
#include "H_Aperture.h"
using namespace std;

H_Aperture::H_Aperture() {
        type = NONE;
        setApertureString();
        x1 = 0;
        x2 = 0;
        x3 = 0;
        x4 = 0;
        fx = 0;
        fy = 0;
}

H_Aperture::H_Aperture(const int dtype, const float size1, const float size2, const float size3, const float size4, const float posx, const float posy) {
    /// @param dtype defines the aperture shape
        /// @param size1, size2, size3, size4 are the geometrical sizes (length/width or great/small radii) in m
	/// @param posx, posy are the (x,y) coordinates of the center of the aperture [m]
        type = dtype;
        setApertureString();
        x1 = size1;
        x2 = size2;
	x3 = size3;
	x4 = size4;
        fx = posx;
        fy = posy;
}

H_Aperture::H_Aperture(const H_Aperture& ap) {
        type = ap.type;
	aptypestring = ap.aptypestring;
        x1 = ap.x1;
        x2 = ap.x2;
        x3 = ap.x3;
        x4 = ap.x4;
        fx = ap.fx;
        fy = ap.fy;
}

H_Aperture& H_Aperture::operator=(const H_Aperture& ap) {
	if(this==&ap) return *this;
	type = ap.type;
	aptypestring = ap.aptypestring;
	x1 = ap.x1;
	x2 = ap.x2;
	x3 = ap.x3;
	x4 = ap.x4;
	fx = ap.fx;
	fy = ap.fy;
	return *this;
}

void H_Aperture::printProperties() const {
        cout << "Aperture shape:" << getTypeString() << ", parameters"<<x1<<", "<<x2<<", "<<x3<<", "<<x4<<endl;
        cout << " \t Center : " << fx << "," << fy << endl;
        return;
}

void H_Aperture::setPosition(const float posx, const float posy) {
	/// @param posx, posy are the (x,y) coordinates of the center of the apertures [m]
	// posx, posy, fx, fy = [m]
	fx = posx;
	fy = posy;
	return;
}

bool H_Aperture::isInside(const float x, const float y) const {
        /// @param x, y are the (x,y) coordinates of the proton, in [m]
        cout << "aperture::isInside" << endl;
        return false;
}


void H_Aperture::setApertureString() {
	switch (type) {
		case NONE: aptypestring = NONENAME; break;
 		case RECTANGULAR: aptypestring = RECTANGULARNAME; break;
 		case ELLIPTIC: aptypestring = ELLIPTICNAME; break;
 		case CIRCULAR: aptypestring = CIRCULARNAME; break;
 		case RECTELLIPSE: aptypestring = RECTELLIPSENAME; break;
		default: aptypestring = NONENAME; break;
	}
}
