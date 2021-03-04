//
#ifndef G__VERTEXFIT_H
#define G__VERTEXFIT_H
//
#include <TMath.h>
#include <TVectorD.h>
#include <TMatrixDSym.h>
#include "ObsTrk.h"
#include <iostream>
//
// Class for vertex fitting

class VertexFit {
	//
	// Vertex fitting with track parameters steering
	// Author: F. Bedeschi, INFN-Pisa, Italy
	// February 10, 2021
	//
private:
	//
	// Inputs
	Int_t fNtr;					// Number of tracks
	TVectorD** fPar;			// Input parameter array
	TMatrixDSym** fCov;			// Input parameter covariances
	// Constraints
	Bool_t fVtxCst;				// Vertex constraint flag
	TVectorD fxCst;				// Constraint value
	TMatrixDSym fCovCst;		// Constraint covariance
	//
	// Results
	Bool_t fVtxDone;			// Flag vertex fit completed
	TVectorD fXv;				// Found vertex
	TMatrixDSym fcovXv;			// Vertex covariance
	Double_t fChi2;				// Vertex fit Chi2
	TVectorD fChi2List;			// List of Chi2 contributions
	//
	// Transient arrays
	Double_t* ffi;				// Fit phases
	TVectorD** fx0i;			// Track expansion points
	TVectorD** fai;				// dx/dphi
	Double_t* fa2i;				// a'Wa
	TMatrixDSym** fDi;			// W-WBW
	TMatrixDSym** fWi;			// (ACA')^-1
	TMatrixDSym** fWinvi;		// ACA'
	//
	// Service routines
	Double_t FastRv1(TVectorD p1, TVectorD p2);		// Fast vertex radius determination
	Double_t FastRv(TVectorD p1, TVectorD p2);		// Fast vertex radius determination
	TMatrixDSym RegInv3(TMatrixDSym& Smat0);		// Regularized 3D matrix inversion 
	TMatrixD Fill_A(TVectorD par, Double_t phi);	// Derivative of track position wrt track parameters
	TVectorD Fill_a(TVectorD par, Double_t phi);	// Derivative of track position wrt track phase
	TVectorD Fill_x0(TVectorD par);					// Track position at dma to z-axis
	TVectorD Fill_x(TVectorD par, Double_t phi);	// Track position at given phase
	void VertexFinder();							// Vertex finder routine
public:
	//
	// Constructors
	VertexFit();										// Initialize waiting for tracks
	VertexFit(Int_t Ntr, ObsTrk** tracks);				// Initialize with ObsTrk tracks
	VertexFit(Int_t Ntr, TVectorD** trkPar, TMatrixDSym** trkCov);	// Initialize with parameters and covariances
	// Destructor
	~VertexFit();
	//
	// Accessors also trigger calculations when needed
	Int_t GetNtrk() { return fNtr; };
	TVectorD GetVtx();
	TMatrixDSym GetVtxCov();
	Double_t GetVtxChi2();
	TVectorD GetVtxChi2List();
	//
	// Handle tracks/constraints
	void AddVtxConstraint(TVectorD xv, TMatrixDSym cov);	// Add gaussian vertex constraint
	void AddTrk(TVectorD par, TMatrixDSym Cov);			// Add track to input list
	void RemoveTrk(Int_t iTrk);							// Remove iTrk track
	//
};

#endif
