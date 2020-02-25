#ifndef G__SOLGRIDCOV_H
#define G__SOLGRIDCOV_H

#include <TVectorD.h>
#include <TMatrixDSym.h>

class SolGeom;

// Class to create geometry for solenoid geometry

class SolGridCov{
  // Class to handle storing and retrieving/interpolation of covariance matrices
private:
  Int_t fNpt;        // Number of pt points in grid
  TVectorD fPta;     // Array of pt points in GeV
  Int_t fNang;       // Number of angle points in grid
  TVectorD fAnga;    // Array of angle points in degrees
  TMatrixDSym *fCov; // Pointers to grid of covariance matrices
  // Service routines
  Int_t GetMinIndex(Double_t xval, Int_t N, TVectorD x); // Find bin
  TMatrixDSym MakePosDef(TMatrixDSym NormMat); // Force positive definitness
public:
  SolGridCov();
  ~SolGridCov();

  void Calc(SolGeom *G);

  // Covariance interpolation
  Double_t GetMinPt()  { return fPta(0); }
  Double_t GetMaxPt()  { return fPta(fNpt - 1); }
  Double_t GetMinAng() { return fAnga(0); }
  Double_t GetMaxAng() { return fAnga(fNang - 1); }
  TMatrixDSym GetCov(Double_t pt, Double_t ang);
};

#endif
