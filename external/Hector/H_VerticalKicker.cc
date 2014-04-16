/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_VerticalKicker.cc
/// \brief Classes aiming at simulating vertical kickers in the beamline

// fk [rad] for kickers !!!!

// local #includes
#include "H_VerticalKicker.h"
#include "H_TransportMatrices.h"

void H_VerticalKicker::setMatrix(const float eloss, const float p_mass, const float p_charge) const {
	extern int kickers_on;
	if(kickers_on) {
		*element_mat = vkickmat(element_length,fk,eloss,p_mass,p_charge);
	} else{
		*element_mat = driftmat(element_length);
	}
	return ;
}
