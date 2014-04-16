#ifndef _H_HorizontalQuadrupole_
#define _H_HorizontalQuadrupole_

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_HorizontalQuadrupole.h
/// \brief Class describing horizontally focussing quadrupoles.

// local #includes
#include "H_Quadrupole.h"

/// Horizontally focussing quadrupoles.
class H_HorizontalQuadrupole : public H_Quadrupole {

	public:
	///     Constructors and destructor
	//@{
		H_HorizontalQuadrupole():H_Quadrupole(HQUADRUPOLE,0.,0.,0.) {init();}
		H_HorizontalQuadrupole(const double s, const double k, const double l) : H_Quadrupole(HQUADRUPOLE,s,k,l){init();}
		H_HorizontalQuadrupole(const string nameE, const double s, const double k, const double l) : H_Quadrupole(nameE,HQUADRUPOLE,s,k,l){init();}
		~H_HorizontalQuadrupole() {return;};
	//@}
	private:
		virtual void setTypeString() {typestring = HQUADRUPOLENAME;};
		virtual void setMatrix(const float, const float, const float) const ;
};

#endif
