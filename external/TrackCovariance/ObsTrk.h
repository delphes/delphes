#ifndef G__OBSTRK_H
#define G__OBSTRK_H

#include <TVector3.h>
#include <TVectorD.h>
#include <TMatrixDSym.h>

class SolGridCov;

// Class to handle smearing of generated charged particle tracks

// Author: F. Bedeschi
//         INFN - Sezione di Pisa, Italy
//
class ObsTrk{
  // Class to handle simulation of tracking resolution
  // Prefix Obs marks variables after resolution smearing
  // Prefix Gen marks variables before resolution smearing
private:
  Double_t fB;      // Solenoid magnetic field
  SolGridCov *fGC;  // Covariance matrix grid
  Double_t fGenQ;   // Generated track charge
  Double_t fObsQ;   // Observed  track charge
  TVector3 fGenX;   // Generated track origin (x,y,z)
  TVector3 fObsX;   // Observed  track origin (x,y,z) @ track min. approach
  TVector3 fGenP;   // Generated track momentum at track origin
  TVector3 fObsP;   // Observed  track momentum @ track minimum approach
  TVectorD fGenPar; // Generated helix track parameters (D, phi0, C, z0, cot(th))
  TVectorD fObsPar; // Observed  helix track parameters (D, phi0, C, z0, cot(th))
  TMatrixDSym fCov; // INterpolated covariance of track parameters
public:
  // x(3) track origin, p(3) track momentum at origin, Q charge, B magnetic field in Tesla
  ObsTrk(TVector3 x, TVector3 p, Double_t Q, Double_t B, SolGridCov *GC); // Initialize and generate smeared track
  ~ObsTrk();
  // Service routines
  TVectorD XPtoPar(TVector3 x, TVector3 p, Double_t Q);
  TVectorD GenToObsPar(TVectorD gPar, SolGridCov *GC);
  TVector3 ParToX(TVectorD Par);
  TVector3 ParToP(TVectorD Par);
  Double_t ParToQ(TVectorD Par);
  // Accessors
  // Generator level X, P, Q
  Double_t GetGenQ() { return fGenQ; }
  TVector3 GetGenX() { return fGenX; }
  TVector3 GetGenP() { return fGenP; }
  // D, phi0, C, z0, cot(th)
  TVectorD GetGenPar() { return fGenPar; }
  // Observed level X, P, Q
  Double_t GetObsQ() { return fObsQ; }
  TVector3 GetObsX() { return fObsX; }
  TVector3 GetObsP() { return fObsP; }
  // D, phi0, C, z0, cot(th)
  TVectorD GetObsPar() { return fObsPar; }
  TMatrixDSym GetCov() { return fCov; }
};

#endif
