#ifndef _H_EllipticAperture_
#define _H_EllipticAperture_

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

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
		H_EllipticAperture(const float, const float, const float, const float);
		~H_EllipticAperture() {return;};
		//@}
		/// Checks whether the point is inside the aperture or not
		virtual bool isInside(const float, const float) const;
		/// Draws the aperture shape.
		virtual void draw() const;
		virtual void printProperties() const;
};

#endif
