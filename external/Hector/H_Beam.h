#ifndef _H_Beam_
#define _H_Beam_

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_Beam.h
/// \brief Describes a set a particles as a beam
///

// c++ #includes
//#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>

// ROOT #includes
////#include "TH2F.h"
////#include "TGraphErrors.h"
////#include "TMultiGraph.h"

// local #includes
#include "H_BeamParticle.h"
#include "H_AbstractBeamLine.h"
#include "H_Parameters.h"
#include "H_OpticalElement.h"
using namespace std;

/// \brief Beam made from a STL vector of H_BeamParticle.
///
/// Example of usage :
/// <PRE>
///   H_Beam mybeam;
///   mybeam.createBeamParticles(100);
///   H_AbstractBeamLine * beamline = new H_AbstractBeamLine(200);
///	... (<-- fills the beam)
///   mybeam.computePath(beamline)
/// </PRE>
class H_Beam {

	public:
	/// Constructors, destructor and operator
	//@{
		H_Beam();
		H_Beam(const H_Beam&);
		H_Beam& operator=(const H_Beam&);
		~H_Beam();
	//@}
	/// Fills the beam with particles
	//@{
		void createBeamParticles(const unsigned int , const double , const double );
		void createBeamParticles(const unsigned int);
	//@}
	/// Fills the beam with particles with incremental offset and no initial transverse momentum
	//@{
		void createXScanningBeamParticles(const unsigned int, const float);
		void createYScanningBeamParticles(const unsigned int, const float);
	//@}
	/// Fills the beam with particles with incremental initial angle and no offset.
	//@{
		void createTXScanningBeamParticles(const unsigned int, const float);
		void createTYScanningBeamParticles(const unsigned int, const float);
	//@}
	/// Returns the ith particle
	//@{
		const H_BeamParticle * getBeamParticle(const unsigned int ) const;
		H_BeamParticle * getBeamParticle(const unsigned int );
	//@}
	/// Adds a particle to the beam
		void add(const H_BeamParticle&);
	/// Compute the position of each particle in the beamline
	//@{	
		void computePath(const H_AbstractBeamLine *, const bool);
		void computePath(const H_AbstractBeamLine *);
	//@}
	// Photon emission by the particle
	// @{
                void emitGamma(const double, const double, const double, const double);
		void emitGamma(const double, const double);
	//@}
	// Propagates the beam until a given s
		void propagate(const float );
	/// Returns the \f$ \beta \f$ function of the beam
	//@{
		float getBetaX(const float , float& );
                float getBetaY(const float , float& );
	//@}
	/// Draws the \f$ \beta \f$ function of the beam
	//@{
        ////	TGraphErrors  * getBetaX(const float, const unsigned int);
	////	TGraphErrors  * getBetaY(const float, const unsigned int);
	//@}
        /// Returns the position of the beam
	//@{
                float getX(const float , float& );
                float getY(const float , float& );
	//@}
	/// Returns the emittance \f$ \epsilon \f$ of the beam in x and y
	//@{
		inline const float getEmittanceX() const {
		    if(!x_disp*tx_disp) cout<<"Warning : Degenerate Beam : x-emittance = 0"<<endl; 	
			return x_disp * tan(tx_disp/URAD)/URAD;
		}
        inline const float getEmittanceY() const { 
			if(!y_disp*ty_disp) cout<<"Warning : Degenerate Beam : y-emittance = 0"<<endl;
			return y_disp * tan(ty_disp/URAD)/URAD;
		}
	//@}
	/// Sets the initial parameters of the beam \f$ s , x , y , \theta_x , \theta_y \f$
	//@{
		void setS(const float fs) { fs_ini = fs;}
		void setX(const float fx) { fx_ini = fx;}
		void setY(const float fy) { fy_ini = fy;}
		void setTX(const float tx) { tx_ini = tx;}
		void setTY(const float ty) { ty_ini = ty;}
		void setE(const float fe) {fe_ini = fe;}
		void setPosition(const float fx, const float fy, const float tx, const float ty, const float fs) {setS(fs); setX(fx); setY(fy); setTX(tx); setTY(ty);}
	//@}
	/// Sets the initial divergence of the beam \f$ \sigma_x , \sigma_y , \sigma (\theta_x), \sigma(\theta_y) \f$
	//@{
		void setDS(const float ds) {s_disp = ds;}
		void setDX(const float dx) {x_disp = dx;}
		void setDY(const float dy) {y_disp = dy;}
		void setDTX(const float dtx) {tx_disp = dtx;}
		void setDTY(const float dty) {ty_disp = dty;}
		void setDE(const float de) {e_disp = de;}
		void setDispersion(const float dx, const float dy, const float dtx, const float dty, const float ds) {setDS(ds); setDX(dx); setDY(dy); setDTX(dtx); setDTY(dty);}
	//@} 
	/// Returns the number of particles which have been stopped
		unsigned int getStoppedNumber(const H_AbstractBeamLine *);
	/// Returns the list of the stopping elements in the beamline
		void getStoppingElements(const H_AbstractBeamLine *, vector<H_OpticalElement>&, vector<int>&);
	/// Prints the initial parameters
		void printInitialState() const;
	/// Prints the properties for each particle
		void printProperties() const;
	/// Prints the list of the stopping elements in the beamline
		void printStoppingElements(const vector<H_OpticalElement>&, const vector<int>&) const;
	/// Returns the number of particle is this beam
		const int getNumberOfBeamParticles() const {return Nparticles;}
	/// Draws the beam profile at a given s
	////	TH2F * drawProfile(const float);
	/// Draws the beam width and height
	//@{
	////	TMultiGraph * drawBeamX(const int) const;
	////	TMultiGraph * drawBeamY(const int) const ;
	//@}

	protected:
	/// List of particles
		vector<H_BeamParticle> beamParticles;
	/// IP position
	//@{
		float fx_ini, fy_ini, fs_ini;
	//@}
	/// IP \f$ \theta \f$ (angular) position
	//@{
		float tx_ini, ty_ini;
	//@}
	/// nominal beam energy
		float fe_ini;
	/// dispersion in position
	//@{
		float x_disp, y_disp, s_disp;
	//@}
	/// dispersion in \f$ \theta \f$
	//@{
		float tx_disp, ty_disp;
	//@}
	/// dispersion in energy
		float e_disp;
	/// Number of particles in this beam
		unsigned int Nparticles;
};

#endif
