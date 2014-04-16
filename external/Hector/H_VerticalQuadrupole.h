#ifndef _H_VerticalQuadrupole_
#define _H_VerticalQuadrupole_

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_VerticalQuadrupole.h
/// \brief Vertically focussing quadrupoles.

// local #includes
#include "H_Quadrupole.h"

/// Vertically focussing quadrupoles.
class H_VerticalQuadrupole : public H_Quadrupole {

	public:
	///     Constructors and destructor
	//@{
		H_VerticalQuadrupole():H_Quadrupole(VQUADRUPOLE,0.,0.,0.) {init();}
		H_VerticalQuadrupole(const double s, const double k, const double l):H_Quadrupole(VQUADRUPOLE,s,k,l) {init();}
		H_VerticalQuadrupole(string nameE, const double s, const double k, const double l):H_Quadrupole(nameE,VQUADRUPOLE,s,k,l) {init();}
		~H_VerticalQuadrupole() {return;};
	//@}
	private:
		virtual void setTypeString() {typestring = VQUADRUPOLENAME;} ;
		virtual void setMatrix(const float, const float, const float) const;
};

#endif
