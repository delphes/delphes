#ifndef _H_RectEllipticAperture_
#define _H_RectEllipticAperture_

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

/// \file H_RectEllipticAperture.h
/// \brief Defines the Rect-Elliptic aperture of beamline elements.

// local #includes
#include "H_Aperture.h"

/// Rect-ellipse apertures
class H_RectEllipticAperture: public H_Aperture {

	public:
		///     Constructors and Destructor
		//@{
		H_RectEllipticAperture():H_Aperture(RECTELLIPSE,0,0,0,0,0,0) {}
		H_RectEllipticAperture(const float,const float,const float,const float, const float, const float);
		~H_RectEllipticAperture() {};
		H_RectEllipticAperture* clone() const;
		//@}
		/// Checks whether the point is inside the aperture or not
		virtual bool isInside(const float, const float) const;
		/// Draws the aperture shape.
		virtual void draw(const float scale=1) const;
	friend std::ostream& operator<< (std::ostream& os, const H_RectEllipticAperture& ap);
};

#endif
