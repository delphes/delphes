#ifndef _H_SectorDipole_
#define _H_SectorDipole_

  /* * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                         *
*                   --<--<--  A fast simulator --<--<--     *
*                 / --<--<--     of particle   --<--<--     *
*  ----HECTOR----<                                          *
*                 \ -->-->-- transport through -->-->--     *
*                   -->-->-- generic beamlines -->-->--     *
*                                                           *
* JINST 2:P09005 (2007)                                     *
*      X Rouby, J de Favereau, K Piotrzkowski (CP3)         *
*       http://www.fynu.ucl.ac.be/hector.html               *
*                                                           *
* Center for Cosmology, Particle Physics and Phenomenology  *
*              Universite catholique de Louvain             *
*                 Louvain-la-Neuve, Belgium                 *
 *                                                         *
   * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/// \file H_SectorDipole.h
/// \brief Classes aiming at simulating sector dipoles.

#include "H_Dipole.h"

/// Sector dipoles.
class H_SectorDipole : public H_Dipole {

	public:
		/// Constructors and destructor
		//@{
		H_SectorDipole():H_Dipole(SDIPOLE,0.,0.,0.) {init();}
		H_SectorDipole(const double s, const double k, const double l) :H_Dipole(SDIPOLE,s,k,l){init();}
		H_SectorDipole(const string& nameE, const double s, const double k, const double l) :H_Dipole(nameE,SDIPOLE,s,k,l){init();}
		~H_SectorDipole() {};
		//@}
		H_SectorDipole* clone() const ;
	private:
		virtual void setTypeString() { typestring = SDIPOLENAME;};
		virtual void setMatrix(const float ,const float, const float) ;
};

#endif
