#ifndef _H_RectangularCollimator_
#define _H_RectangularCollimator_

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

/// \file H_RectangularCollimator.h
/// \brief Class describing Rectangular Collimators

// local #includes
#include "H_Drift.h"

/// R-Collimators for the LHC beamline.
class H_RectangularCollimator : public H_Drift {

	public:
	/// Constructors and destructor
	//@{
		H_RectangularCollimator():H_Drift() {type = RCOLLIMATOR; init();}
		H_RectangularCollimator(const double, const double );
		H_RectangularCollimator(const string&, const double, const double );
		~H_RectangularCollimator() {};
	//@}
		H_RectangularCollimator* clone() const;
		void init();
	private:
	 	virtual void setTypeString() {typestring=RCOLLIMATORNAME;};
		virtual void setMatrix(const float, const float, const float) ;
	
	friend std::ostream& operator<< (std::ostream& os, const H_RectangularCollimator& el);
};

#endif
