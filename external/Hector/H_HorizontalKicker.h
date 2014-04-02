#ifndef _H_HorizontalKicker_
#define _H_HorizontalKicker_

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

/// \file H_HorizontalKicker.h
/// \brief Classes aiming at simulating horizontal kickers in beamline.
///
/// fk [rad] for kickers !!!!

// local #includes
#include "H_Kicker.h"

/// Horizontal kickers for the beamline.
class H_HorizontalKicker : public H_Kicker {

	public:
	/// Constructors and destructor
	//@{
		H_HorizontalKicker():H_Kicker(HKICKER,0.,0.,0.) {init();}
		H_HorizontalKicker(const double s, const double k, const double l) :H_Kicker(HKICKER,s,k,l){init();}
		H_HorizontalKicker(const string &nameE, const double s, const double k, const double l) :H_Kicker(nameE,HKICKER,s,k,l){init();}
		~H_HorizontalKicker() {};
	//@}
		H_HorizontalKicker* clone() const;
	private:
		virtual void setTypeString() {typestring=HKICKERNAME;};
		virtual void setMatrix(const float, const float, const float) ;
};

#endif
