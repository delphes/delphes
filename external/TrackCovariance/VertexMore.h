//
#ifndef G__VERTEXMORE_H
#define G__VERTEXMORE_H
//
#include <TMath.h>
#include <TVectorD.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TMatrixDSym.h>
#include "TrkUtil.h"
#include "VertexFit.h"
#include <vector>
#include <iostream>
//
// Class to include found vertices in vertex fitting
//
//	Author: F. Bedeschi, INFN-Pisa, Italy
//

class VertexMore: public TrkUtil
{
	//
	// Vertex fitting with vertices
	// Author: F. Bedeschi, INFN-Pisa, Italy
	// April 7, 2022
	//
private:
	//
	// Inputs
	VertexFit* fV;						// Input vertex
	TVectorD fXv;						// running vertex position
	TMatrixDSym fXvCov;					// running vertex error matrix
	Bool_t fUnits;						// Default is m, can switch to mm
	Double_t fa;						// conversion constant C = a/(2Pt)
	Int_t fNtr;						// NUmber of tracks in vertex
	TVectorD fPar;						// Vertex track parameters
	TMatrixDSym fCov;					// Vertex track covariance
	TVectorD MakeVpar();					// Init vertex parameters
	TMatrixDSym MakeVcov();					// Init Vertex parameter covariance
	std::vector<Double_t> fQ;				// Track Charges
	std::vector<TVector3*> fpi;				// Track Momenta
	std::vector<TMatrixDSym*> fCpi;				// Track Momentum errors
	TVector3 fP;						// Total momentum
	Double_t fQtot;						// Vertex charge
	TMatrixDSym fCp;					// Total momentum errors
	TVectorD fBigPar;					// Full vertex/momenta vector
	TMatrixDSym fBigCov;					// Full covariance matrix of vertex position and momenta
	void Mprt(TMatrixD M, Bool_t Opt);			// Compact matrix printing
	void Init(VertexFit* V);				// Constructor initializations
	void FillBigCov();					// Fill fBigCov
	void FillBigPar();					// Fill fBigPar
	//	
	Double_t dSdD(Int_t i);					// Derivative of phase wrt D
	Double_t dSdC(Int_t i);					// Derivative of phase wrt C
	void CalcParCov();					// Calculate parameters and covariance
	TMatrixD dPdAlf(Int_t i);			// Derivatives of momentum wrt parameters
	TVectorD dPds(Int_t i);				// Derivative of momentum wrt phase
	TMatrixD dPdX(Int_t i);				// Derivatives of momentum wrt vertex position
	TMatrixD dXdAlf(Int_t i);			// Derivatives of position wrt parameters
	TMatrixD DparDx(TVector3 xv, TVector3 pv, Double_t Q); // Derivatives of parameters wrt point on track
	TMatrixD DparDp(TVector3 xv, TVector3 pv, Double_t Q); // Derivatives of parameters wrt track momentum
	TMatrixD DpDa0(Int_t i, Int_t k);		// Derivatives of momentum wrt initial track parameters
	TMatrixD BuildPPcov(Int_t i, Int_t j);		// Build <dPi,dPj> covariance from VertexFit info
	TMatrixD BuildPXcov(Int_t i);			// Build <dPi,dXv> covariance from VertexFit info
	//
	// Mass constraints
	//
	Int_t fNc;					// Nr. of mass constraints
	TVectorD fMassC;				// Constraint vector
	TMatrixD fMderC;				// Constraint derivatives wrt momenta	
	std::vector<Double_t> fMass;			// Mass of ith constraint
	std::vector<Int_t>    fNtracks;			// nr. tracks in ith constraint
	std::vector<Double_t*> fmasses;			// Input masses lists
	std::vector<Int_t*>    flists;			// Input list of tracks
	void UpdateConstr(TVectorD p);			// Update constraints and derivatives
public:
	//
	// Constructors
	VertexMore(VertexFit *V);			// Initialize with found vertex (use meter as units)
	VertexMore(VertexFit* V, Bool_t Opt);		// Initialize with unit option TRUE = use mm, FALSE = use meters
	// Destructor
	~VertexMore();
	//
	VertexFit* GetVinput(){ return fV; };		// Return pointer to initial vertex
	TVectorD GetVpar(){ return fPar;};		// Get vertex track parameters
	TMatrixDSym GetVcov(){ return fCov;};		// Get vertex track covariance
	Double_t GetCharge(Int_t i) { return fQ[i]; };
	TVector3 GetMomentum(Int_t i) { return *fpi[i]; };		// Momentum of track i at vertex
	TMatrixDSym GetMomentumC(Int_t i) { return *fCpi[i]; };		// Momentum errors of track i at vertex
	TVector3 GetTotalP() { return fP; };				// Total vertex momentum
	Double_t GetTotalQ() { return fQtot; };				// Total vertex charge
	TMatrixDSym GetTotalPcov() { return fCp; };			// Total vertex momentum errors
	TMatrixDSym GetBigCov() { return fBigCov; };			// Vertex/momenta full covariance
	TVectorD GetBigPar(){ return fBigPar;};				// Vertex/momenta array
	TVectorD GetXv(){ return fXv;};					// Vertex position
	TMatrixDSym GetXvCov(){ return fXvCov;};			// Vertex covariance;
	//
	// Mass constraints
	void AddMassConstraint(Double_t Mass, Int_t Ntr, Double_t* masses, Int_t* list); 
	void MassConstrFit();
};

#endif
