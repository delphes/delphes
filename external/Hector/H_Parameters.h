#ifndef _Hector_parameters_
#define _Hector_parameters_

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_Parameters.h
/// \brief Class aiming at gathering all parameters that must be defined.
///
/// Units : angles [\f$ \mu \f$rad], distances [\f$ \mu \f$m], energies [GeV], c=[1].

/* from physics and maths */
   /// proton mass [GeV]
const double MP=0.93827;
   /// proton charge [e]
const double QP=1;
   /// pi
#define PI 3.14159265359
   /// conversion factor for \f$\mu\f$rad <-> rad
#define URAD 1000000.
   /// in thx = thx* / 3 <==> thx* = THX thx
#define THX 3.86 
   /// in dy = 10 thy*   <==> thy* = dy / THY
#define THY 51.4       

/* beam parameters */
   /// beam nominal energy in GeV
//#define BE  7000.
const double BE=7000.;
   /// beam energy divergence, in GeV
#define SBE 0.79
   /// beam nominal energy in TeV
#define BETEV 7.
   /// beam S @ IP
#define PS 0.
   /// beam X @ IP
#define PX  -500.      
   /// beam Y @ IP
#define PY  0.      
   /// beam longitudinal dispersion
#define SS 0.
   /// beam lateral width SX @ IP
#define SX  16.63
// #define SX  0.      
   /// beam lateral width SY @ IP
#define SY  16.63
// #define SY 0.
   /// beam transverse direction angle TX @ IP
#define TX 0.
   /// beam transverse direction angle TY @ IP
#define TY 0.      
   /// beam angular divergence STX @ IP
//#define STX 0.
#define STX 30.23      
   /// beam angular divergence STY @ IP
//#define STY 0.
#define STY 30.23      
   /// beam dispersion
//#define D   120000.
const double D=120000.;
	/// half crossing angle at IP [\f$ \mu \f$RAD]
#define CRANG 142.5 

/* roman pots parameters */
   /// granularity in X position
#define GRANPOSX 5.
   /// granularity in Y position
#define GRANPOSY 5.
   /// granularity in X angle
#define GRANANGX 10.
   /// granularity in Y angle
#define GRANANGY 10.
   /// Distance between rp's
#define DISTRP 4000000.
   /// RP resolution in X, for smearing
#define RESX 10.
   /// RP resolution in Y, for smearing
#define RESY 10.
   /// Radius of the hole in the RP
#define RADIUS 1000.

/* display parameters */
   /// Verbose mode ?
#define VERBOSE 0 

/// include Pythia libraries ? (not included on some ROOT installations) 
//#define _include_pythia_

#endif

