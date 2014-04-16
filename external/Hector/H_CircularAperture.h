#ifndef _H_CircularAperture_
#define _H_CircularAperture_

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_CircularAperture.h
/// \brief Defines the circular aperture of beamline elements.

// local #includes
#include "H_EllipticAperture.h"


/// Circular apertures
class H_CircularAperture: public H_EllipticAperture {

	public:
		/// Constructors and destructor
		//@{
		H_CircularAperture():H_EllipticAperture(0,0,0,0) {type = CIRCULAR; setApertureString();}
		H_CircularAperture(const float, const float, const float);
		~H_CircularAperture() {return;};
		//@}
		virtual void printProperties() const;
};

#endif
