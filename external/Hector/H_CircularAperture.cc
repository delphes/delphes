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

/// \file H_CircularAperture.cc
/// \brief Defines the circular aperture of beamline elements.

// C++ #includes
#include <iostream>

// C #includes
#include <cmath> // needed for fabs and log

// local #includes
#include "H_CircularAperture.h"
using namespace std;

H_CircularAperture* H_CircularAperture::clone() const {
	return new H_CircularAperture(x1,fx,fy);
}

std::ostream& operator<< (std::ostream& os, const H_CircularAperture& ap) {
        os     << "Aperture shape:" << ap.aptypestring << ", aperture radius : " << ap.x1 << endl;
        os     << " \t Center : " << ap.fx << "," << ap.fy << endl;
        return os;
}
