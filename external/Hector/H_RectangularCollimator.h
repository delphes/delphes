#ifndef _H_RectangularCollimator_
#define _H_RectangularCollimator_

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_RectangularCollimator.h
/// \brief Class describing Rectangular Collimators

// local #includes
#include "H_OpticalElement.h"

/// R-Collimators for the LHC beamline.
class H_RectangularCollimator : public H_OpticalElement {

	public:
	/// Constructors and destructor
	//@{
		H_RectangularCollimator():H_OpticalElement(RCOLLIMATOR,0.,0.,0.) {init();}
		H_RectangularCollimator(const double, const double );
		H_RectangularCollimator(const string, const double, const double );
		~H_RectangularCollimator() { return; };
	//@}
		virtual void printProperties() const;
		void init();

	private:
	 	virtual void setTypeString() {typestring=RCOLLIMATORNAME;};
		virtual void setMatrix(const float, const float, const float) const ;

};

#endif
