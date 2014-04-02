#ifndef _H_Marker_
#define _H_Marker_

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

/// \file H_Marker.h
/// \brief Class defining a marker in the beamline (e.g. for interaction point)

// local includes
#include "H_Drift.h"

/// \brief Class defining a marker in the beamline (e.g. for interaction point)
class H_Marker : public H_Drift {

	public:
	/// Constructors and destructor
	//@{
		H_Marker():H_Drift() { type = MARKER; init();}
        	H_Marker(const double s):H_Drift(s,0.) { type =MARKER; init();}
	        H_Marker(const string& nameE, const double s):H_Drift(nameE,s,0.) { type=MARKER; init();}
	        ~H_Marker() { };
	//@}
		H_Marker* clone() const ;
	        void init();

	private:
        	virtual void setTypeString() {typestring = MARKERNAME;};
	        virtual void setMatrix(const float , const float, const float) ;
	friend std::ostream& operator<< (std::ostream& os, const H_Marker& el);
};

#endif
