#ifndef _H_Drift_
#define _H_Drift_

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_Drift.h
/// \brief Class aiming at simulating LHC beam drift.

// local includes
#include "H_OpticalElement.h"

/// \brief Class aiming at simulating LHC beam drift. 
class H_Drift : public H_OpticalElement {

	public:
	/// Constructors and destructor
	//@{
		H_Drift():H_OpticalElement(DRIFT,0.,0.,0.) {init();}
        	H_Drift(const double s, const double l):H_OpticalElement(DRIFT,s,0.,l){init();}
	        H_Drift(const string nameE, const double s, const double l):H_OpticalElement(nameE,DRIFT,s,0.,l){init();}
	        ~H_Drift() { return; };
	//@}
        	virtual void printProperties() const;
	        void init();

	private:
        	virtual void setTypeString() {typestring = DRIFTNAME;};
	        virtual void setMatrix(const float, const float, const float) const ;
};

#endif
