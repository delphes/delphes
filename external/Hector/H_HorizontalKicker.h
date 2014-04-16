/// \file H_HorizontalKicker.h
/// \brief Classes aiming at simulating horizontal kickers in beamline.
///
/// fk [rad] for kickers !!!!

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

#ifndef _H_HorizontalKicker_
#define _H_HorizontalKicker_

// local #includes
#include "H_Kicker.h"

/// Horizontal kickers for the beamline.
class H_HorizontalKicker : public H_Kicker {

	public:
	/// Constructors and destructor
	//@{
		H_HorizontalKicker():H_Kicker(HKICKER,0.,0.,0.) {init();}
		H_HorizontalKicker(const double s, const double k, const double l) :H_Kicker(HKICKER,s,k,l){init();}
		H_HorizontalKicker(const string nameE, const double s, const double k, const double l) :H_Kicker(nameE,HKICKER,s,k,l){init();}
		~H_HorizontalKicker() {return;};
	//@}
	private:
		virtual void setTypeString() {typestring=HKICKERNAME;};
		virtual void setMatrix(const float, const float, const float) const ;
};

#endif
