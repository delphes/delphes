#ifndef _H_RectangularAperture_
#define _H_RectangularAperture_

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

/// \file H_RectangularAperture.h
/// \brief Defines the rectangular aperture of beamline elements.

// local #includes
#include "H_Aperture.h"

/// Rectangular apertures
class H_RectangularAperture : public H_Aperture {

	public:
		/// Constructors and destructor
		//@{
		H_RectangularAperture():H_Aperture(RECTANGULAR,0,0,0,0,0,0) {}
		H_RectangularAperture(const float,const float,const float,const float);
		~H_RectangularAperture() {};
		H_RectangularAperture* clone() const;
		//@}
		/// Checks whether the point is inside the aperture or not 
		virtual bool isInside(const float, const float) const;
		/// Draws the aperture shape.
		virtual void draw(const float scale=1) const;
	friend std::ostream& operator<< (std::ostream& os, const H_RectangularAperture& ap);
};

#endif
