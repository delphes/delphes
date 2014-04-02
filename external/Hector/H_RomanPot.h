#ifndef _H_RomanPot_
#define _H_RomanPot_

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

/// \file H_RomanPot.h
/// \brief Roman pot class

// local #includes
#include "H_Drift.h"

#define RP_LENGTH 0.0001
#define RP_HEIGHT 10000
/// Roman pot as an optics element.
class H_RomanPot : public H_Drift {

	public:
	///     Constructors and destructor
	//@{
		H_RomanPot():H_Drift() {type = RP; init();}
		H_RomanPot(const string&, const double, const double);
		H_RomanPot(const double, const double);
		~H_RomanPot() {};
	//@}
		H_RomanPot* clone() const ;
		void init();
	private:
		virtual void setTypeString() {typestring=RPNAME;};
		virtual void setMatrix(const float, const float, const float) ;
	friend std::ostream& operator<< (std::ostream& os, const H_RomanPot& el);
};

#endif
