#ifndef _H_BeamParticle_
#define _H_BeamParticle_

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_BeamParticle.h
/// \brief Class aiming at simulating a particle in the beamline

// from IP to RP, with emission of a photon of defined energy and Q.
// Units : angles [rad], distances [m], energies [GeV], masses [GeV], c=[1].
// !!! no comment statement at the end of a #define line !!!

// c++ #includes
#include <vector>

// ROOT #includes
#include "TMatrixD.h"
#include "TVectorD.h"
//#include "TPolyLine3D.h"
#include "TRandom.h"

// local #includes
#include "H_Parameters.h"
#include "H_AbstractBeamLine.h"
#include "H_OpticalElement.h"

using namespace std;

// local defines
#define LENGTH_VEC 5
#define INDEX_X 0
#define INDEX_TX 1
#define INDEX_Y 2
#define INDEX_TY 3
#define INDEX_S 4
// (x,theta_x,y,theta_y,s)

/// Defines a particle from the beam and its transport through the beamline
class H_BeamParticle {

	public:
		void init();
		/// Constructors and Destructor
		//@{
		H_BeamParticle();
		H_BeamParticle(const H_BeamParticle&);
		H_BeamParticle(const double, const double);
		H_BeamParticle& operator=(const H_BeamParticle&);
		~H_BeamParticle() {delete stop_position; if(!stop_element) delete stop_element; positions.clear(); return; }
		//@}
		/// Smears the (x,y) coordinates of the particle [\f$ \mu m \f$]
		void smearPos(const double dx=SX,const double dy=SY, TRandom* r=gRandom);
		/// Smears the (x,y) angular coordinates of the particle [\f$ \mu rad \f$]
		void smearAng(const double tx=STX, const double ty=STY, TRandom* r=gRandom);
		/// Smears the Energy of the particle [GeV]
		void smearE(const double erre=SBE, TRandom* r=gRandom);
		/// Smears the longitudinal position of the particle [\f$ \mu m \f$]
		void smearS(const double errs=SS, TRandom* r=gRandom);
		/// Sets the energy [GeV].
		void setE(const double);
		/// Sets the particle 4-momentum \f$ P^\mu \f$
		void set4Momentum(const double, const double, const double, const double);	
		/// Clears H_BeamParticle::positions and sets the initial one.
		void setPosition(const double , const double , const double , const double , const double );
		/// Returns the particle mass [GeV]
		double getM() const {return mp;};
		/// Returns the particle charge [e]
		double getQ() const {return qp;};
		/// Returns the current x coordinate [\f$ \mu \f$m]
		double getX() const {return fx;};
		/// Returns the current y coordinate [\f$ \mu \f$m]
		double getY() const {return fy;};
		/// Returns the current s coordinate [m]
		inline double getS() const {return fs;};
		/// Returns the current \f$ \theta_x \f$ angular coordinate [\f$ \mu \f$rad]
		inline double getTX() const {return thx;};
		/// Returns the current \f$ \theta_y \f$ angular coordinate [\f$ \mu \f$rad]
		inline double getTY() const {return thy;};
		/// Returns the current particle energy [GeV]
		inline double getE() const {return energy;};
		/// Returns all the positions 
		vector<TVectorD> getPositions() const {return positions;};
		bool isPhysical() const {return isphysical;};
		/// \brief Simulates the emission of a photon in a random direction
		///
	        /// For \f$ p_{1} \rightarrow p_{2} \gamma \f$, kinematics imposes that
	        /// \f$ Q^{2} = E^{2}_{\gamma} -p^{2}_{1} -p^{2}_{2} + 2p_{1}p_{2} cos(\theta)  \f$ where \f$ \theta \f$ is the particle scattering angle and \f$ p_{i} = \|\vec{p_{i}}\| \f$. <BR>
        	/// So, \f$ Q^{2}_{min} = E^{2}_{\gamma} - (p_{1}+p_{2})^{2} \f$ and \f$ Q^{2}_{max} = E^{2}_{\gamma} - (p_{1}-p_{2})^{2} \f$.<BR>
	        /// As \f$ E^{2}_{\gamma} - (p_{1}-p_{2})^{2} \f$ could be numerically instable, we use here another form of this formula : <BR>
        	/// \f$ Q^{2}_{max} = -2 * \big( \frac{M_{p} E_{\gamma}}{p_{1}+p_{2}} \big) \big[ 1 + \frac{E^{2}_{1} + E^{2}_{2} - M^{2}_{p} }{ E_{1} E_{2} + p_{1} p_{2}} \big]    \f$
		//@{
		void emitGamma(const double, const double, const double, const double);
		void emitGamma(const double, const double);
		//@}
                /// uses Pythia to generate some inelastic pp->pX collision as background
		void doInelastic();
		/// \brief Propagates the particule across the beamline until the s coordinate is reached
		///
		/// Caution : "computePath" should be used before any "propagate" call <BR>
		/// Caution : "stopped" is not included in "propagate" : please run it afterward if needed 
		void propagate(const double ) ;
		/// Propagates the particle accross the beamline until a given element
		void propagate(const H_AbstractBeamLine *, const H_OpticalElement *);
		/// Propagates the particle accross the beamline until a given element
		void propagate(const H_AbstractBeamLine *, const string);
		/// Propagates the particle until the end of the beamline
		void propagate(const H_AbstractBeamLine *);
		/// Returns the phase vector of the particle
		const TMatrixD * getV() const;
		/// Returns the current phase vector of the particle (in H_BeamParticle::positions)
		const TVectorD * getPosition(const int ) const;
		/// Prints the properties of the particle
		void printProperties() const;
		/// Prints the phase vector of the particle
		void printV() const;
		/// Returns the element where the particle has been stopped
		const H_OpticalElement * getStoppingElement() const;
		/// Checks if the particle has been stopped in any element of the beamline
		bool stopped(const H_AbstractBeamLine *);
		/// Returns the StopPosition vector
		inline const TVectorD * getStopPosition() const { return stop_position; };
		// returns (-1,-1,-1,-1,-1) if not stopped (and then hasstopped is false)
		/// Shows all the vectors \f$ (x, \theta_x, y, \theta_y ,s) \f$ in H_BeamParticle::positions
		void showPositions() const; 
		/// Returns the particle path in the beamline
		////TGraph * getPath(const int , const int ) const;
		/// Draws the particle path in the beamline in 3D
		////TPolyLine3D * getPath3D(const H_AbstractBeamLine *, const bool, const int, const int) const;
		/// Computes the position of the particle at the end of each element of the beam, without non linear effects
		void computePath(const H_AbstractBeamLine *);
		/// Computes the position of the particle at the end of each element of the beam.
		void computePath(const H_AbstractBeamLine *, const bool);
		/// Computes the position of the particle at the end of each element of the beam.
		void computePath(const H_AbstractBeamLine &, const bool);
		/// Clears H_BeamParticle::positions but keeps the initial vector.
		void resetPath();

	private:
		/// Particle mass [GeV]
		double mp;
		/// Particle charge [e]
		double qp;
		/// Longitudinal co-moving coordinate [m]
		double fs;
		/// Transverse (horizontal) coordinate [m]
		double fx;
		/// Transverse (vertical) coordinate [m]
		double fy;
		/// Direction of the 3-momentum in the horizontal plane [rad]
		double thx;
		/// Direction of the 3-momentum in the vertical plane [rad]
		double thy;
		/// Kinetic energy of the particle [GeV]
		double energy;
		/// True if the particle has stopped (i.e. : if the particle transverse position has been out of any optics element aperture). <BR> See H_BeamParticle::stopped method. <BR>Default = false.
		bool hasstopped;
		/// True if the particle has lost some (E,Q), i.e. if H_BeamParticle::emitGamma was used. Default = false.
		bool hasemitted;
		/// False if the particle has emitted a photon with impossible (E,Q). Default = true.
		bool isphysical;
		/// Vector (x,tx,y,ty,s) where the particle has stopped.
		TVectorD * stop_position; 
		/// Element of the beamline (H_OpticalElement) where the particle has stopped. 
		H_OpticalElement * stop_element;
		/// List of (x,tx,y,ty,s) vectors, after each optical element of the beam
		vector<TVectorD> positions; // vector (x,tx,y,ty,s) after each optical element of the beam ([m],[rad],[m],[rad],[m])
		/// Adds a new vector (x,tx,y,ty,s) at the end of H_BeamParticle::positions
  		void addPosition(const double , const double , const double , const double , const double );
};
#endif
