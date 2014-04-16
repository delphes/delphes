#ifndef _H_RectangularAperture_
#define _H_RectangularAperture_

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

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
		~H_RectangularAperture() {return;};
		//@}
		/// Checks whether the point is inside the aperture or not 
		virtual bool isInside(const float, const float) const;
		/// Draws the aperture shape.
		virtual void draw() const;
		virtual void printProperties() const;
};

#endif
