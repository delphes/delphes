#ifndef _H_Quadrupole_
#define _H_Quadrupole_

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_Quadrupole.h
/// \brief Class describing quadrupoles.

// local #includes
#include "H_OpticalElement.h"

/// Abstract class aiming at simulating LHC beam quadrupoles.
class H_Quadrupole : public H_OpticalElement {

	public:
	///     Constructors and destructor
	//@{
		H_Quadrupole():H_OpticalElement() {}
		H_Quadrupole(const int dtype, const double s, const double k, const double l) : H_OpticalElement(dtype,s,k,l) {}
		H_Quadrupole(const string nameE, const int dtype, const double s, const double k, const double l) : H_OpticalElement(nameE,dtype,s,k,l) {}
		virtual ~H_Quadrupole() {return;};
	//@}	
		virtual void printProperties() const;
		void init();

 	protected:
		virtual void setTypeString() =0;
		virtual void setMatrix(const float, const float, const float) const =0;
	   	
};

#endif
