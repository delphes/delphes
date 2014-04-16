#ifndef _H_VerticalKicker_
#define _H_VerticalKicker_

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_VerticalKicker.h
/// \brief Classes aiming at simulating vertical kickers in the beamline.
///
/// fk [rad] for kickers !!!!


// local #includes
#include "H_Kicker.h"

/// Vertical kickers for the beamline.
class H_VerticalKicker : public H_Kicker {

	public:
	/// Constructors and destructor
	//@{
		H_VerticalKicker():H_Kicker(VKICKER,0.,0.,0.) {init();}
		H_VerticalKicker(const double s, const double k, const double l) :H_Kicker(VKICKER,s,k,l){init();}
		H_VerticalKicker(const string nameE, const double s, const double k, const double l) :H_Kicker(nameE,VKICKER,s,k,l){init();}
		~H_VerticalKicker() {return;};
	//@}
	private:
		virtual void setTypeString() {typestring=VKICKERNAME;};
		virtual void setMatrix(const float, const float, const float) const ;
};

#endif
