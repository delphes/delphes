#ifndef _H_BeamLine_
#define _H_BeamLine_

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
		H_BeamLine();
		H_BeamLine(const H_BeamLine& );
		H_BeamLine(const int, const float);
		H_BeamLine& operator=(const H_BeamLine& );
		~H_BeamLine() {};
	//@
		///     Finds the IP position (s) from the MAD table. Should be "IP5" or "IP1".
		void findIP(const string& filename, const string& ipname="IP5");
		///     Reader for the external MAD table
		void fill(const string& filename, const int dir=1, const string& ipname="IP5");
		///     Returns the IP position (s)
		double getIP() const {return ips;};
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

