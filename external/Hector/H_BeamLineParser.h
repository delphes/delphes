#ifndef _H_BeamLineParser_
#define _H_BeamLineParser_

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_BeamLineParser.h
/// \brief Reader for madx tables
///
/// Notes : Verifier que tous les SBEND sont toujours appeles MB. Et seulement comme ca !
/// Verifier qu'il n'y a pas de probleme d'inversion H/V QUADRUPOLES
/// no distinction between H and Vkickers ?
/// The identification of the element is based on the values of k1l, k2l, hkick, vkick and on their name.

// c++ #includes
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
using namespace std;

// local defines
#define MADX_UNKNOWN 0
#define MADX_NAME    1
#define MADX_KEYWORD 2
#define MADX_S       3
#define MADX_L       4
#define MADX_K0L     5
#define MADX_K1L     6
#define MADX_K2L     7
#define MADX_K3L     8
#define MADX_HKICK   9
#define MADX_VKICK  10
#define MADX_BETX   11
#define MADX_ALFX   12
#define MADX_MUX    13
#define MADX_DX     14
#define MADX_DPX    15
#define MADX_X      16
#define MADX_PX     17
#define MADX_BETY   18
#define MADX_ALFY   19
#define MADX_MUY    20
#define MADX_DY     21
#define MADX_DPY    22
#define MADX_Y      23
#define MADX_PY     24
#define MADX_APERTYPE 25
#define MADX_APER_1 26
#define MADX_APER_2 27
#define MADX_APER_3 28
#define MADX_APER_4 29
#define MADX_KICK 30
#define MADX_PARENT 31
/* Apertures # parameters
	CIRCLE 1
	ELLIPSE	2
	RECTANGLE	2
	LHCSCREEN	3
	MARGUERITE	3
	RECTELLIPSE	4
	RACETRACK	3
*/

extern int column_identification(const string );

/// \brief Reader for madx tables to use in H_BeamLine
///
/// More info on http://project-mad9.web.cern.ch/project-mad9/mad/mad9/user/index.html
class H_BeamLineParser {

	public:
	/// Constructor and destructor
	//@{
		H_BeamLineParser() {init();}
		~H_BeamLineParser() {return;}
	//@}
		void init();
		/// Retrieve the data from the line it reads, and sets the corresponding variable.
		void setProperties(istream& , const unsigned int );
		void printProperties() const;
		string name, apertype, keyword, parent;
		/// Optical element longitudinal (co-moving) coordinate.
		double s;
		/// Length of the element.
		double l;
	/// Magnetic field strength.
	//@{
		double k0l, k1l, k2l, k3l, hkick, vkick;
	//@}
	/// Phase (mu), \f$ \alpha \f$(alf) and \f$ \beta \f$(bet) functions, given by MAD.
	//@{
		double mux, muy, betx, alfx, bety, alfy;
	//@}
	/// Positions and their dispersion
	//@{ 
		double x, y, dx, dy;
	//@}
	/// Momentum and their dispersion
	//@{
		double px, py, dpx, dpy;
	//@}
	/// Aperture parameters
	//@{
		double aper_1, aper_2, aper_3, aper_4;
	//@}
};

#endif
