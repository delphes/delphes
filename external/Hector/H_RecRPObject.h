#ifndef _H_RecRPObject_
#define _H_RecRPObject_

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_RecRPObject.h
/// \brief Reconstructed information from objects detected by two roman pots 

// local #includes
#include "H_Parameters.h"
#include "H_AbstractBeamLine.h"
#include <math.h>
using namespace std;

#define NOT_YET_COMPUTED -666
// trivial (TM), angle compensation (AM) and position compensation (PM) methods
#define TM 1
#define AM 2
#define PM 3

/// Reconstructed information from objects detected by two roman pots
class H_RecRPObject {
/// Uses (s1,x1,y1) and (s2,x2,y2) to compute the energy, the virtuality of the event, as well as the positions (x,y)  and angles (thx, thy) at IP.
	public:
		/// Constructors, destructor and operators
		//@{
		H_RecRPObject();
		H_RecRPObject(const float, const float, const H_AbstractBeamLine& );
		H_RecRPObject(const H_RecRPObject&);
		H_RecRPObject& operator=(const H_RecRPObject&);
		~H_RecRPObject() {delete matrp1; delete matrp2; delete thebeam; return;};
		//@}

		/// Getters 
		//@{
		inline float getX1() const {return x1; /*horizontal position at first roman pot, in \mu m */}
		inline float getY1() const {return y1; /*vertical position at first roman pot, in \mu m */}
		inline float getS1() const {return s1; /*longitudinal position of the first roman pot, in m, from IP*/}
		inline float getX2() const {return x2; /*horizontal position at second RP in \mu m*/}
		inline float getY2() const {return y2; /*vertical position at second RP in \mu m*/}
		inline float getS2() const {return s2; /*longitudinal position of the second RP, in m, from IP*/}
		inline float getTXRP() const {return URAD*atan((x2-x1)/((s2-s1)*URAD)); /* horizontal angle at first RP, in \mu rad*/ }
		inline float getTYRP() const {return URAD*atan((y2-y1)/((s2-s1)*URAD)); /* vertical angle at first RP, in \mu rad*/ }
		float getX0() {return (x0==NOT_YET_COMPUTED) ? computeX0():x0; /*reconstructed horizontal position at IP*/}
		float getY0() {return (y0==NOT_YET_COMPUTED) ? computeY0():y0; /*reconstructed vertical position at IP*/}
		float getTXIP() {return (thx==NOT_YET_COMPUTED) ? computeTX():thx; /*reconstructed horizontal angle at IP*/}
		float getTYIP() {return (thy==NOT_YET_COMPUTED) ? computeTY():thy; /*reconstructed vertical angle at IP*/}
		float getE(const unsigned int); /*returns the reconstructed energy*/
		float getE(); /*returns the reconstructed energy if already computed*/ 
		float getQ2() {return (virtuality==NOT_YET_COMPUTED) ? computeQ2():virtuality; /*returns the reconstructed virtuality*/}
		//@}

		// Sets the proton hit positions
		void setPositions(const float, const float, const float, const float);
		// Shows the variable content.
		void printProperties() const;		
	protected:
		/// Measured particle coordinates at RP (X - horizontal and Y - vertical in [\f$ \mu \f$m], S -longitudinal in [m])
		//@{
		float x1, x2, y1, y2, s1, s2; 
		//@}

		/// Reconstructed positions and angles at IP in [\f$ \mu m\f$] and [\f$ \mu rad\f$] 
		//@{
		float x0, y0, thx, thy;
		//@}
	
		/// Reconstructed energy and virtuality at IP in GeV and GeV\f$^2\f$ 
		//@{
		float energy, virtuality;
		//@}
	
		float computeX0();
		float computeY0();	
		float computeTX();
		float computeTY();
		/// Energy reconstruction : trivial method
		float computeE_TM();
		/// Energy reconstruction : angle compensation method
		float computeE_AM();
		/// Energy reconstruction : position compensation method
		float computeE_PM();
		/// Virtuality reconstruction. Energy should be reconstructed before.
		float computeQ2();
		/// Calibrates the energy reconstruction with respect to the chromaticity of the transfer matrices
		float getECorrectionFactor(const unsigned int, const unsigned int );
		/// The beamline : 
		H_AbstractBeamLine * thebeam;
		/// The matrices
		//@{
		TMatrix * matrp1;
		TMatrix * matrp2;
		//@}
		/// The correction factors
		//@{
		float corr1_TM, corr2_TM, corr1_AM, corr2_AM;
		//@}
};

#endif
