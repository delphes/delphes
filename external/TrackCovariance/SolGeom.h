#ifndef G__SOLGEOM_H
#define G__SOLGEOM_H

#include <Rtypes.h>
#include <TString.h>

// Class to create geometry for solenoid geometry
// Simplified implementations with cylindrical and disk layers
//
// Author: F. Bedeschi, INFN - Pisa

class SolGeom{
  //
  // Units are m
  //
private:
  const Int_t fNlMax = 700; // Maximum number of layers

  // B field
  Double_t fB; // B field in Tesla

  // Barrel layer properties
  Int_t fNlay;       // Total number of layers
  Int_t fBlay;       // Number of barrel layers
  Int_t fFlay;       // Number of forward/backward layers
  Int_t fNm;         // Nr. measurement layers
  Int_t *ftyLay;     // Layer type 1 = R (barrel) or 2 = z (forward/backward)
  TString *fLyLabl;  // Layer label
  // Barrel: PIPE, VTXLOW, VTXHIGH, DCHCANI, DCH, DCHCANO, BSILWRP, MAG, BPRESH
  // Fw/Bw: VTXDSK, DCHWALL, FSILWRP, FRAD, FPRESH
  Double_t *fxMin;   // Minimum dimension z for barrel  or R for forward
  Double_t *fxMax;   // Maximum dimension z for barrel  or R for forward
  Double_t *frPos;   // R/z location of layer
  Double_t *fthLay;  // Thickness (meters)
  Double_t *frlLay;  // Radiation length (meters)
  Int_t    *fnmLay;  // Number of measurements in layers (1D or 2D)
  Double_t *fstLayU; // Stereo angle (rad) - 0(pi/2) = axial(z) layer - Upper side
  Double_t *fstLayL; // Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
  Double_t *fsgLayU; // Resolution Upper side (meters) - 0 = no measurement
  Double_t *fsgLayL; // Resolution Lower side (meters) - 0 = no measurement
  Bool_t   *fflLay;  // measurement flag = T, scattering only = F
//
//
  Double_t fRmin;	// Radius of first barrel layer
  Double_t fZminPos;	// Z of first disk in positive direction
  Double_t fZminNeg;	// Z of first disk in negative direction
  void SetMinBoundaries();	// define inner box for fast simulation

public:
  SolGeom();
  ~SolGeom();

  void Read(const char *data);
  void SetBz(const Double_t Bz);
  
  Double_t B()                     { return fB; }
  Int_t    Nl()                    { return fNlay; }
  Int_t    Nm()                    { return fNm; }
  Int_t    NBl()                   { return fBlay; }
  TString  lLabl(Int_t nlayer)     { return fLyLabl[nlayer]; }
  Int_t    lTyp(Int_t nlayer)      { return ftyLay[nlayer]; }
  Double_t lxMin(Int_t nlayer)     { return fxMin[nlayer]; }
  Double_t lxMax(Int_t nlayer)     { return fxMax[nlayer]; }
  Double_t lPos(Int_t nlayer)      { return frPos[nlayer]; }
  Double_t lTh(Int_t nlayer)       { return fthLay[nlayer]; }
  Double_t lX0(Int_t nlayer)       { return frlLay[nlayer]; }
  Int_t    lND(Int_t nlayer)       { return fnmLay[nlayer]; }
  Double_t lStU(Int_t nlayer)      { return fstLayU[nlayer]; }
  Double_t lStL(Int_t nlayer)      { return fstLayL[nlayer]; }
  Double_t lSgU(Int_t nlayer)      { return fsgLayU[nlayer]; }
  Double_t lSgL(Int_t nlayer)      { return fsgLayL[nlayer]; }
  Bool_t   isMeasure(Int_t nlayer) { return fflLay[nlayer]; }
  //
  // Define cylindrical box to use for fast simulation
  //
  Double_t GetRmin() { return fRmin; }
  Double_t GetZminPos() { return fZminPos; }
  Double_t GetZminNeg() { return fZminNeg; }
};

#endif
