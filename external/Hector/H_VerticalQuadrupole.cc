/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_VerticalQuadrupole.cc
/// \brief Vertically focussing quadrupoles

// local #includes
#include "H_VerticalQuadrupole.h"

void H_VerticalQuadrupole::setMatrix(const float eloss, const float p_mass, const float p_charge) const {
	if (fk<0)  { if(VERBOSE) cout<<"\t ERROR : k1 should be > 0 for H_VerticalQuadrupole (" << name << ")!"<<endl; }
		if (fk !=0 ) *element_mat = vquadmat(element_length,fk,eloss, p_mass,p_charge);
		else  {
			*element_mat = driftmat(element_length);
			if(VERBOSE) cout<<"\t WARNING : k1= 0, drift-like quadrupole (" << name << ") !" << endl;
   		}
	return ;
}
