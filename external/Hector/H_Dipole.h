/// \file H_Dipole.h
/// \brief Class aiming at simulating LHC beam dipoles.

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

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
		H_Dipole(const string nameE, const int dtype, const double s, const double k, const double l):H_OpticalElement(nameE,dtype,s,k,l){}
	        virtual ~H_Dipole() {return;};
		//@}
		/// Prints the properties of the element
        	virtual void printProperties() const;
	        void init();

 	protected:
    		virtual void setTypeString() =0;
	    	virtual void setMatrix(const float, const float, const float) const =0;
    	
};

#endif
