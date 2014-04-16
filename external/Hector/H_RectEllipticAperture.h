#ifndef _H_RectEllipticAperture_
#define _H_RectEllipticAperture_

/// \file H_RectEllipticAperture.h
/// \brief Defines the Rect-Elliptic aperture of beamline elements.

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

// local #includes
#include "H_Aperture.h"

/// Rect-ellipse apertures
class H_RectEllipticAperture: public H_Aperture {

	public:
		///     Constructors and Destructor
		//@{
		H_RectEllipticAperture():H_Aperture(RECTELLIPSE,0,0,0,0,0,0) {}
		H_RectEllipticAperture(const float,const float,const float,const float, const float, const float);
		~H_RectEllipticAperture() {return;};
		//@}
		virtual void printProperties() const;
		/// Checks whether the point is inside the aperture or not
		virtual bool isInside(const float, const float) const;
		/// Draws the aperture shape.
		virtual void draw() const;
};

#endif
