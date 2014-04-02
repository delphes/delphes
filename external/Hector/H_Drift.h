#ifndef _H_Drift_
#define _H_Drift_

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
	        H_Drift(const string& nameE, const double s, const double l):H_OpticalElement(nameE,DRIFT,s,0.,l){init();}
	        ~H_Drift() { };
	//@}
	        void init();
		virtual void printProperties() const { cout << *this; return;};
		H_Drift* clone() const;

	protected:
        	virtual void setTypeString() {typestring = DRIFTNAME;};
	        virtual void setMatrix(const float, const float, const float) ;
	friend std::ostream& operator<< (std::ostream& os, const H_Drift& el);
};

#endif
