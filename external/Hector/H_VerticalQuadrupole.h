#ifndef _H_VerticalQuadrupole_
#define _H_VerticalQuadrupole_

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
		H_VerticalQuadrupole(const string& nameE, const double s, const double k, const double l):H_Quadrupole(nameE,VQUADRUPOLE,s,k,l) {init();}
	~H_VerticalQuadrupole() {};
//@}
	H_VerticalQuadrupole* clone() const ;
private:
		virtual void setTypeString() {typestring = VQUADRUPOLENAME;} ;
		virtual void setMatrix(const float, const float, const float) ;
};

#endif
