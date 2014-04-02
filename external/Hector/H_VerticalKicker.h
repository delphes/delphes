#ifndef _H_VerticalKicker_
#define _H_VerticalKicker_

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
		H_VerticalKicker(const string& nameE, const double s, const double k, const double l) :H_Kicker(nameE,VKICKER,s,k,l){init();}
		~H_VerticalKicker() {};
	//@}
		H_VerticalKicker* clone() const ;
	private:
		virtual void setTypeString() {typestring=VKICKERNAME;};
		virtual void setMatrix(const float, const float, const float) ;
};

#endif
