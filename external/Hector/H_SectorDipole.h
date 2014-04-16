/// \file H_SectorDipole.h
/// \brief Classes aiming at simulating sector dipoles.

#ifndef _H_SectorDipole_
#define _H_SectorDipole_

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

#include "H_Dipole.h"

/// Sector dipoles.
class H_SectorDipole : public H_Dipole {

	public:
		/// Constructors and destructor
		//@{
		H_SectorDipole():H_Dipole(SDIPOLE,0.,0.,0.) {init();}
		H_SectorDipole(const double s, const double k, const double l) :H_Dipole(SDIPOLE,s,k,l){init();}
		H_SectorDipole(const string nameE, const double s, const double k, const double l) :H_Dipole(nameE,SDIPOLE,s,k,l){init();}
		~H_SectorDipole() {return;};
		//@}
	private:
		virtual void setTypeString() { typestring = SDIPOLENAME;};
		virtual void setMatrix(const float ,const float, const float) const ;
};

#endif
