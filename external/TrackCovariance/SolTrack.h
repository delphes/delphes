//
#ifndef G__SOLTRK_H
#define G__SOLTRK_H
//
#include <TMath.h>
#include <TVector3.h>
#include <TMatrixDSym.h>
#include "SolGeom.h"
#include "TrkUtil.h"
#include <TGraph.h>
//
//
// Class to store track information
// Assumes that the geometry has been initialized
//
class SolTrack: public TrkUtil
{
	// 
	// Track handling class
	// Assume tracks originate from (0,0) for the time being
	//
private:
	Int_t fNl;	// Actual number of layers
	SolGeom *fG;	// Geometry
	Double_t fp[3];	// px, py, pz momentum 
	Double_t fx[3];	//  x,  y,  z track origin 
	Double_t fpar[5];	// D, phi0, C, z0, cot(theta)
	//
	TMatrixDSym fCov;		// Full covariance matrix
	//
	//
public:
	//
	// Constructors
	SolTrack(Double_t *x, Double_t *p, SolGeom *G);
	SolTrack(TVector3 x, TVector3 p, SolGeom* G);
	SolTrack(Double_t D, Double_t phi0, Double_t C, Double_t z0, Double_t ct, SolGeom *G);
	// Destructor
	~SolTrack();
	// Accessors
	// Position (at minimum approach)
	Double_t x() { return fx[0]; }
	Double_t y() { return fx[1]; }
	Double_t z() { return fx[2]; }
	// Momentum  (at minimum approach)
	Double_t px() { return fp[0]; }
	Double_t py() { return fp[1]; }
	Double_t pz() { return fp[2]; }
	Double_t pt() { return TMath::Sqrt(fp[0] * fp[0] + fp[1] * fp[1]); }
	Double_t p()  { return TMath::Sqrt(fp[0] * fp[0] + fp[1] * fp[1] + fp[2] * fp[2]); }
	// Track parameters
	Double_t D()    { return fpar[0]; }
	Double_t phi0() { return fpar[1]; }
	Double_t C()    { return fpar[2]; }
	Double_t z0()   { return fpar[3]; }
	Double_t ct()   { return fpar[4]; }
	// Covariance
	TMatrixDSym Cov()		{ return fCov; }
	// Track parameter covariance calculation
	void CovCalc(Bool_t Res, Bool_t MS);
	// Parameter errors
	Double_t s_D()    { return TMath::Sqrt(fCov(0, 0)); }
	Double_t s_phi0() { return TMath::Sqrt(fCov(1, 1)); }
	Double_t s_C()    { return TMath::Sqrt(fCov(2, 2)); }
	Double_t s_pt()   { return 2 * s_C()*pt() / (0.2998*fG->B()); }	// Dpt/pt
	Double_t s_z0()   { return TMath::Sqrt(fCov(3, 3)); }
	Double_t s_ct()   { return TMath::Sqrt(fCov(4, 4)); }
	//
	// Track hit management
	Int_t nHit();
	Int_t nmHit();
	Bool_t HitLayer(Int_t Layer, Double_t &R, Double_t &phi, Double_t &zz);
	Int_t HitList(Int_t *&ihh, Double_t *&rhh, Double_t *&zhh);
	Int_t HitListXYZ(Int_t *&ihh, Double_t *&Xh, Double_t *&Yh, Double_t *&Zh);
	Int_t FirstHit(Double_t &Xfirst, Double_t &Yfirst, Double_t &Zfirst);	// First hit position
	//
	// Track graph
	TGraph *TrkPlot();	// Graph with R-z plot of track trajectory
	//
	// Make normalized matrix positive definite
	TMatrixDSym MakePosDef(TMatrixDSym NormMat);
};
//
#endif
