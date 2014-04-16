/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_SectorDipole.cc
/// \brief Classes aiming at simulating sector dipoles

#include "H_SectorDipole.h"
#include "H_TransportMatrices.h"

void H_SectorDipole::setMatrix(const float eloss, const float p_mass, const float p_charge) const {
	if (fk !=0 ) *element_mat = sdipmat(element_length,fk,eloss,p_mass,p_charge);
	else  {
		*element_mat = driftmat(element_length);
		if(VERBOSE) cout<<"\t WARNING : k0= 0, drift-like dipole (" << name << ") !" << endl;
   		}
	return ;
}
