#ifndef _Hector_parameters_
#define _Hector_parameters_

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

/// \file H_Parameters.h
/// \brief Class aiming at gathering all parameters that must be defined.
///
/// Units : angles [\f$ \mu \f$rad], distances [\f$ \mu \f$m], energies [GeV], c=[1].

#include <cmath>

/* from physics and maths */
   /// proton mass [GeV]
const double MP=0.93827;
   /// proton charge [e]
const double QP=1.;
   /// conversion factor for \f$\mu\f$rad <-> rad
const double URAD=1000000.;

/* beam parameters */
   /// beam nominal energy in GeV
const double BE=7000.;
   /// beam energy divergence, in GeV
const double SBE=0.79;
   /// beam nominal energy in TeV
const double BETEV=7.;
   /// beam S @ IP
const double PS=0.;
   /// beam X @ IP
const double PX=-500.;    
   /// beam Y @ IP
const double PY=0.;
   /// beam longitudinal dispersion
const double SS=0.;
   /// beam lateral width SX @ IP
const double SX=16.63;
   /// beam lateral width SY @ IP
const double SY=16.63;
   /// beam transverse direction angle TX @ IP
const double TX=0.;
   /// beam transverse direction angle TY @ IP
const double TY=0.;
   /// beam angular divergence STX @ IP
const double STX=30.23;
   /// beam angular divergence STY @ IP
const double STY=30.23;
	/// half crossing angle at IP [\f$ \mu \f$RAD]
const double CRANG=142.5;

// local defines, used in H_BeamParticle & H_OpticalElements
enum {INDEX_X=0, INDEX_TX, INDEX_Y, INDEX_TY, INDEX_S, LENGTH_VEC};
// (x,theta_x,y,theta_y,s)

/// include Pythia libraries ? (not included on some ROOT installations)
//#define _include_pythia_

const unsigned int TM = 0; // not used anymore. left for backward compatibility
const unsigned int AM = 1; // not used anymore. left for backward compatibility


/* display parameters */
   /// Verbose mode ?
const bool VERBOSE=false; 

#endif

