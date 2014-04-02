#ifndef _H_Quadrupole_
#define _H_Quadrupole_

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

/// \file H_Quadrupole.h
/// \brief Class describing quadrupoles.

// local #includes
#include "H_OpticalElement.h"

/// Abstract class aiming at simulating LHC beam quadrupoles.
class H_Quadrupole : public H_OpticalElement {

	public:
	///     Constructors and destructor
	//@{
		H_Quadrupole():H_OpticalElement() {}
		H_Quadrupole(const int dtype, const double s, const double k, const double l) : H_OpticalElement(dtype,s,k,l) {}
		H_Quadrupole(const string& nameE, const int dtype, const double s, const double k, const double l) : H_OpticalElement(nameE,dtype,s,k,l) {}
		virtual ~H_Quadrupole() {};
	//@}	
		virtual H_Quadrupole* clone() const =0 ;
		virtual void printProperties() const { cout << *this; return;};
		void init();

 	protected:
		virtual void setTypeString() =0;
		virtual void setMatrix(const float, const float, const float) = 0;
	   	
	friend std::ostream& operator<< (std::ostream& os, const H_Quadrupole& el);
};

#endif
