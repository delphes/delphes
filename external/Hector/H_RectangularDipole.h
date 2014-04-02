#ifndef _H_RectangularDipole_
#define _H_RectangularDipole_

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

/// \file H_RectangularDipole.h
/// \brief Classes aiming at simulating LHC beam rectangular dipoles.

#include "H_Dipole.h"

/// Rectangle dipoles.
class H_RectangularDipole : public H_Dipole {

	public:
		/// constructor
		//@{
		H_RectangularDipole():H_Dipole(RDIPOLE,0.,0.,0.) {init();}
		H_RectangularDipole(const double s, const double k, const double l) :H_Dipole(RDIPOLE,s,k,l){init();}
		H_RectangularDipole(const string& nameE, const double s, const double k, const double l) :H_Dipole(nameE,RDIPOLE,s,k,l){init();}
		~H_RectangularDipole() {};
		//@}
		H_RectangularDipole* clone() const ;
	private:
		virtual void setTypeString() { typestring = RDIPOLENAME;};
		virtual void setMatrix(const float, const float, const float) ;
};

#endif
