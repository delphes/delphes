#ifndef _H_Marker_
#define _H_Marker_

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_Marker.h
/// \brief Class defining a marker in the beamline (e.g. for interaction point)

// local includes
#include "H_OpticalElement.h"

/// \brief Class defining a marker in the beamline (e.g. for interaction point)
class H_Marker : public H_OpticalElement {

	public:
	/// Constructors and destructor
	//@{
		H_Marker():H_OpticalElement(MARKER,0.,0.,0.) {init();}
        	H_Marker(const double s):H_OpticalElement(MARKER,s,0.,0.){init();}
	        H_Marker(const string nameE, const double s):H_OpticalElement(nameE,MARKER,s,0.,0.){init();}
	        ~H_Marker() { return; };
	//@}
        	virtual void printProperties() const;
	        void init();

	private:
        	virtual void setTypeString() {typestring = MARKERNAME;};
	        virtual void setMatrix(const float , const float, const float) const ;
};

#endif
