#ifndef _H_HorizontalQuadrupole_
#define _H_HorizontalQuadrupole_

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
		H_HorizontalQuadrupole(const string& nameE, const double s, const double k, const double l) : H_Quadrupole(nameE,HQUADRUPOLE,s,k,l){init();}
		~H_HorizontalQuadrupole() {};
	//@}
		H_HorizontalQuadrupole* clone() const ;
	private:
		virtual void setTypeString() {typestring = HQUADRUPOLENAME;};
		virtual void setMatrix(const float, const float, const float) ;
};

#endif
