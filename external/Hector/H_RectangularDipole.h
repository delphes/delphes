#ifndef _H_RectangularDipole_
#define _H_RectangularDipole_

/// \file H_RectangularDipole.h
/// \brief Classes aiming at simulating LHC beam rectangular dipoles.

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

#include "H_Dipole.h"

/// Rectangle dipoles.
class H_RectangularDipole : public H_Dipole {

	public:
		/// constructor
		//@{
		H_RectangularDipole():H_Dipole(RDIPOLE,0.,0.,0.) {init();}
        	H_RectangularDipole(const double s, const double k, const double l) :H_Dipole(RDIPOLE,s,k,l){init();}
	        H_RectangularDipole(const string nameE, const double s, const double k, const double l) :H_Dipole(nameE,RDIPOLE,s,k,l){init();}
        	~H_RectangularDipole() {return;};
		//@}
	private:
	        virtual void setTypeString() { typestring = RDIPOLENAME;};
        	virtual void setMatrix(const float, const float, const float) const ;
};

#endif
