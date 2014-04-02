#ifndef _H_EllipticAperture_
#define _H_EllipticAperture_

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

/// \file H_EllipticAperture.h
/// \brief Defines the elliptic aperture of beamline elements.

// local #includes
#include "H_Aperture.h"

/// Elliptic apertures
class H_EllipticAperture: public H_Aperture {

	public:
		/// Constructors and destructors
		//@{
		H_EllipticAperture():H_Aperture(ELLIPTIC,0,0,0,0,0,0) {}
		H_EllipticAperture(const int, const float, const float, const float, const float);
		H_EllipticAperture(const float, const float, const float, const float);
		~H_EllipticAperture() {};
		virtual H_EllipticAperture* clone() const;
		//@}
		/// Checks whether the point is inside the aperture or not
		virtual bool isInside(const float, const float) const;
		/// Draws the aperture shape.
		virtual void draw(const float scale=1) const;
	friend std::ostream& operator<< (std::ostream& os, const H_EllipticAperture& ap);
};

#endif
