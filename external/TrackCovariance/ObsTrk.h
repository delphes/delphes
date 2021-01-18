//
#ifndef G__OBSTRK_H
#define G__OBSTRK_H
#include <TVector3.h>
#include <TVectorD.h>
#include <TMatrixDSym.h>
#include <TDecompChol.h>
#include "SolGridCov.h"
//
// Class to handle smearing of generated charged particle tracks
//
// Author: F. Bedeschi
//         INFN - Sezione di Pisa, Italy
//
class ObsTrk{
	//
	// Class to handle simulation of tracking resolution
	// Prefix Obs marks variables after resolution smearing
	// Prefix Gen marks variables before resolution smearing
	//
private:
	Double_t fB;						// Solenoid magnetic field
	SolGridCov *fGC;					// Covariance matrix grid
	Double_t fGenQ;					// Generated track charge
	Double_t fObsQ;					// Observed  track charge
	TVector3 fGenX;					// Generated track origin (x,y,z)
	TVector3 fObsX;					// Observed  track origin (x,y,z) @ track min. approach
	TVector3 fGenP;					// Generated track momentum at track origin
	TVector3 fObsP;					// Observed  track momentum @ track minimum approach
	TVectorD fGenPar;				// Generated helix track parameters (D, phi0, C, z0, cot(th))
	TVectorD fGenParACTS;			// Generated helix track parameters (D, z0, phi0, th, q/p, time
	TVectorD fGenParILC;				// Generated helix track parameters (w, phi0, d0, z0, tan(lambda))
	TVectorD fObsPar;				// Observed  helix track parameters (D, phi0, C, z0, cot(th))
	TVectorD fObsParACTS;			// Observed  helix track parameters (D, z0, phi0, th, q/p, time
	TVectorD fObsParILC;				// Observed  helix track parameters (d0, phi0, w, z0, tan(lambda))
	TMatrixDSym fCov;				// INterpolated covariance of track parameters
	TMatrixDSym fCovACTS;			// Covariance of track parameters in ACTS format
									// (D, z0, phi0, theta, q/p, time)
	TMatrixDSym fCovILC;				// Covariance of track parameters in ILC format
									// (d0, phi0, w, z0, tan(lambda))
	//
	// Conversion to ACTS parametrization
	//
	TVectorD ParToACTS(TVectorD Par);		// Parameter conversion
	TMatrixDSym CovToACTS(TMatrixDSym Cov);	// Covariance
	//
	// Conversion to ILC parametrization
	//
	TVectorD ParToILC(TVectorD Par);		// Parameter conversion
	TMatrixDSym CovToILC(TMatrixDSym Cov);	// Covariance conversion
	//
public:
	//
	// Constructors
	// x(3) track origin, p(3) track momentum at origin, Q charge, B magnetic field in Tesla
	ObsTrk(TVector3 x, TVector3 p, Double_t Q, Double_t B, SolGridCov *GC);	// Initialize and generate smeared track
	ObsTrk(Double_t *x, Double_t *p, Double_t Q, Double_t B, SolGridCov* GC);	// Initialize and generate smeared track
  // Destructor
	~ObsTrk();
	//
	// Service routines
	//
	TVectorD XPtoPar(TVector3 x, TVector3 p, Double_t Q);
	TVectorD GenToObsPar(TVectorD gPar, SolGridCov *GC);
	TVector3 ParToX(TVectorD Par);
	TVector3 ParToP(TVectorD Par);
	Double_t ParToQ(TVectorD Par);
	//
	// Accessors
	//
	// Generator level:
	// X, P, Q
	Double_t GetGenQ()	{ return fGenQ; }
	TVector3 GetGenX()	{ return fGenX; }
	TVector3 GetGenP()	{ return fGenP; }
	// D, phi0, C, z0, cot(th)
	TVectorD GetGenPar()	{ return fGenPar; }
	// D, z0, phi0, theta, q/p, time
	TVectorD GetGenParACTS()	{ return fGenParACTS; }
	// d0, phi0, w, z0, tan(lambda)
	TVectorD GetGenParILC()	{ return fGenParILC; }
	// Observed level X, P, Q
	Double_t GetObsQ()	{ return fObsQ; }
	TVector3 GetObsX()	{ return fObsX; }
	TVector3 GetObsP()	{ return fObsP; }
	// D, phi0, C, z0, cot(th)
	TVectorD GetObsPar()	{ return fObsPar; }
	// D, z0, phi0, theta, q/p, time
	TVectorD GetObsParACTS()	{ return fObsParACTS; }
	// d0, phi0, w, z0, tan(lambda)
	TVectorD GetObsParILC()	{ return fObsParILC; }
	// Covariances
	TMatrixDSym GetCov(){ return fCov; }
	TMatrixDSym GetCovACTS(){ return fCovACTS; }
	TMatrixDSym GetCovILC(){ return fCovILC; }
};

#endif
