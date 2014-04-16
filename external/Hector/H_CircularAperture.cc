/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_CircularAperture.cc
/// \brief Defines the circular aperture of beamline elements.

// C++ #includes
#include <iostream>

// C #includes
#include <cmath> // needed for fabs and log

// local #includes
#include "H_CircularAperture.h"
using namespace std;

/// Circular apertures
H_CircularAperture::H_CircularAperture(const float r, const float posx, const float posy) :H_EllipticAperture(r,r,posx,posy) { 
        /// @param r is the radius of the circular shape
        /// @param posx, posy are the (x,y) coordinates of the center of the circle
	type= CIRCULAR; 
}

void H_CircularAperture::printProperties() const {
        cout << "Aperture shape:" << getTypeString() << ", aperture radius : " << x1 << endl;
        cout << " \t Center : " << fx << "," << fy << endl;
        return;
}
