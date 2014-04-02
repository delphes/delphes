#ifndef _H_CircularAperture_
#define _H_CircularAperture_

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

/// \file H_CircularAperture.h
/// \brief Defines the circular aperture of beamline elements.

// local #includes
#include "H_EllipticAperture.h"


/// Circular apertures
class H_CircularAperture: public H_EllipticAperture {

	public:
		/// Constructors and destructor
		//@{
		H_CircularAperture():H_EllipticAperture(CIRCULAR,0,0,0,0) {};
		H_CircularAperture(const float r, const float posx, const float posy) : H_EllipticAperture(CIRCULAR,r,r,posx,posy) {};
		/// @param r is the radius of the circular shape
        	/// @param posx, posy are the (x,y) coordinates of the center of the circle
		~H_CircularAperture() {};
		H_CircularAperture* clone() const;
		//@}
	friend std::ostream& operator<< (std::ostream& os, const H_CircularAperture& ap);
};

#endif
