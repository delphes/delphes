#ifndef _H_Kicker_
#define _H_Kicker_

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

/// \file H_Kicker.h
/// \brief Classes aiming at simulating kickers in LHC beamline.
/// fk [rad] for kickers !!!!

// local #includes
#include "H_OpticalElement.h"

/// Abstract class aiming at simulating kickers in LHC beamline.
class H_Kicker : public H_OpticalElement {

public:
	/// Constructors and destructor
	//@{
		H_Kicker():H_OpticalElement() {}
		H_Kicker(const int dtype, const double s, const double k, const double l):H_OpticalElement(dtype,s,k,l){}
		H_Kicker(const string& nameE, const int dtype, const double s, const double k, const double l):H_OpticalElement(nameE,dtype,s,k,l){}
		virtual ~H_Kicker() {};
	//@}
		virtual H_Kicker* clone() const = 0;
		virtual void printProperties() const { cout << *this; return;};
		void init();

 	protected:
		virtual void setTypeString() = 0;
		virtual void setMatrix(const float, const float, const float) = 0;
	friend std::ostream& operator<< (std::ostream& os, const H_Kicker& el);
};

#endif
