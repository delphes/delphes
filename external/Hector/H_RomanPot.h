#ifndef _H_RomanPot_
#define _H_RomanPot_

/// \file H_RomanPot.h
/// \brief Roman pot class

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

// local #includes
#include "H_OpticalElement.h"

#define RP_LENGTH 0.0001
#define RP_HEIGHT 10000
/// Roman pot as an optics element.
class H_RomanPot : public H_OpticalElement {

	public:
	///     Constructors and destructor
	//@{
		H_RomanPot():H_OpticalElement(RP,0.,0.,RP_LENGTH) {init();}
		H_RomanPot(const string, const double, const double);
		H_RomanPot(const double, const double);
		~H_RomanPot() { return; };
	//@}
		virtual void printProperties() const;
		void init();

	private:
		virtual void setTypeString() {typestring=RPNAME;};
		virtual void setMatrix(const float, const float, const float) const;

};

#endif
