//
#ifndef G__VERTEXFIT_H
#define G__VERTEXFIT_H
//
#include <TMath.h>
#include <TVectorD.h>
#include <TMatrixDSym.h>
#include "TrkUtil.h"
#include "ObsTrk.h"
#include <vector>
#include <iostream>
//
// Class for vertex fitting

class VertexFit: public TrkUtil
{
	//
	// Vertex fitting with track parameters steering
	// Author: F. Bedeschi, INFN-Pisa, Italy
	// February 10, 2021
	//
private:
	//
	// Inputs
	Int_t fNtr;							// Number of tracks
	std::vector<TVectorD*> fPar;		// Input parameter array
	std::vector<TVectorD*> fParNew;		// Updated parameter array
	std::vector<TMatrixDSym*> fCov;		// Input parameter covariances
	std::vector<TMatrixDSym*> fCovNew;	// Updated parameter covariances
	// Constraints
	Bool_t fVtxCst;					// Vertex constraint flag
	TVectorD fxCst;					// Constraint value
	TMatrixDSym fCovCst;			// Constraint 
	TMatrixDSym fCovCstInv;			// Inverse of constraint covariance
	//
	// Results
	Bool_t fVtxDone;			// Flag vertex fit completed
	Double_t fRstart;			// Starting value of vertex radius (0 = none)
	TVectorD fXv;				// Found vertex
	TMatrixDSym fcovXv;			// Vertex covariance
	Double_t fChi2;				// Vertex fit Chi2
	TVectorD fChi2List;			// List of Chi2 contributions
	//
	// Work arrays
	std::vector<Double_t> ffi;				// Fit phases
	std::vector<TVectorD*> fx0i;			// Track expansion points
	std::vector<TVectorD*> fai;				// dx/dphi
	std::vector<TVectorD*> fdi;				// x-shift
	std::vector<Double_t> fa2i;				// a'Wa
	std::vector<TMatrixD*> fAti;			// A transposed
	std::vector<TMatrixDSym*> fDi;			// W-WBW
	std::vector<TMatrixDSym*> fWi;			// (ACA')^-1
	std::vector<TMatrixDSym*> fWinvi;		// ACA'
	//
	// Service routines
	void ResetWrkArrays();							// Clear work arrays
	//Double_t StartRadius();							// Starting vertex radius determination
	//Double_t FastRv(TVectorD p1, TVectorD p2);		// Fast vertex radius determination
	//TMatrixDSym RegInv(TMatrixDSym& Smat0);		// Regularized 3D matrix inversion 
	//TMatrixD Fill_A(TVectorD par, Double_t phi);	// Derivative of track position wrt track parameters
	//TVectorD Fill_a(TVectorD par, Double_t phi);	// Derivative of track position wrt track phase
	TVectorD Fill_x0(TVectorD par);					// Track position at dma to z-axis
	TVectorD Fill_x(TVectorD par, Double_t phi);	// Track position at given phase
	void UpdateTrkArrays(Int_t i);					// Fill track realted arrays
	void VtxFitNoSteer();							// Vertex fitter routine w/o parameter steering
	void VertexFitter();							// Vertex fitter routine w/  parameter steering
public:
	//
	// Constructors
	VertexFit();										// Initialize waiting for tracks
	VertexFit(Int_t Ntr, ObsTrk** tracks);			// Initialize with ObsTrk tracks
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
	void AddTrk(TVectorD *par, TMatrixDSym *Cov);		// Add track to input list
	void RemoveTrk(Int_t iTrk);				// Remove iTrk track
	void SetStartR(Double_t R) { fRstart = R; };		// Set starting radius
	//
};

#endif
