#ifndef _H_BeamLine_
#define _H_BeamLine_

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_BeamLine.h
/// \brief Class aiming at retrieving features of the real beam optical elements.

// local #includes
#include "H_AbstractBeamLine.h"

// C++ #includes
#include <string>
using namespace std;

// Caution : mad conventions for vertically(horizontally) focusing quadrupoles
// are opposite to Wille's. See fill() for more details.

/// Reads external files and retrieves features of the real beam optical elements.
class H_BeamLine : public H_AbstractBeamLine {

	public:
	///     Constructors and destructor
	//@{
		H_BeamLine():H_AbstractBeamLine() {direction=1; ips=0; };
		H_BeamLine(const H_BeamLine& );
		H_BeamLine(const int, const float);
		H_BeamLine& operator=(const H_BeamLine& );
		~H_BeamLine() {return;};
	//@
		///     Finds the IP position (s) from the MAD table. Should be "IP5" or "IP1".
		//@{
		void findIP(const string);
		void findIP(const string, const string);
		//@}
		///     Reader for the external MAD table
		//@{
		void fill(const string);
		void fill(const string, const int, const string );
		//@}
		///     Returns the IP position (s)
		double getIP() {return ips;};
		///		Returns positions and angles of beam at IP
		double* getIPProperties();

	private:
		int direction; // to or from the IP.
		double ips; // s-position of the IP [m]
		double ipx; // x-position of the IP [µm]
		double ipy; // y-position of the IP [µm]
		double iptx; // x-angle of momentum at IP [µrad]
		double ipty; // y-angle of momentum at IP [µrad]
};

#endif

