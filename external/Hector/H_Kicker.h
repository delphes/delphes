/// \file H_Kicker.h
/// \brief Classes aiming at simulating kickers in LHC beamline.
/// fk [rad] for kickers !!!!

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

#ifndef _H_Kicker_
#define _H_Kicker_

// local #includes
#include "H_OpticalElement.h"

/// Abstract class aiming at simulating kickers in LHC beamline.
class H_Kicker : public H_OpticalElement {

public:
	/// Constructors and destructor
	//@{
		H_Kicker():H_OpticalElement() {}
		H_Kicker(const int dtype, const double s, const double k, const double l):H_OpticalElement(dtype,s,k,l){}
		H_Kicker(const string nameE, const int dtype, const double s, const double k, const double l):H_OpticalElement(nameE,dtype,s,k,l){}
		virtual ~H_Kicker() {return;};
	//@}
	/// prints the kicker properties
		virtual void printProperties() const;
		void init();

 	protected:
		virtual void setTypeString() =0;
		virtual void setMatrix(const float, const float, const float) const =0;
};

#endif
