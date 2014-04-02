/// \file H_Dipole.h
/// \brief Class aiming at simulating LHC beam dipoles.

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

#ifndef _H_Dipole_
#define _H_Dipole_

#include "H_OpticalElement.h"

/// Abstract class aiming at simulating dipoles.
class H_Dipole : public H_OpticalElement {

	public:
		/// Constructors and destructor
		//@{
		H_Dipole():H_OpticalElement() {}
		H_Dipole(const int dtype, const double s, const double k, const double l):H_OpticalElement(dtype,s,k,l){}
		H_Dipole(const string& nameE, const int dtype, const double s, const double k, const double l):H_OpticalElement(nameE,dtype,s,k,l){}
	        virtual ~H_Dipole() {};
		//@}
		/// Prints the properties of the element
		virtual void printProperties() const { cout << *this; return;};
		virtual H_Dipole* clone() const =0;
		void init();

 	protected:
    		virtual void setTypeString() =0;
	    	virtual void setMatrix(const float, const float, const float) = 0;
    
	friend std::ostream& operator<< (std::ostream& os, const H_Dipole& el);
};

#endif
