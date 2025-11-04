//
#ifndef G__KALMANCK_H
#define G__KALMANCK_H
//
#include <iostream>
#include <vector>
#include <TGraph.h>
#include <TAxis.h>
#include <TMath.h>
#include <TVector3.h>
#include <TVector.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TString.h>
#include "examples/classes/SolGeom.h"
#include "external/TrackCovariance/TrkUtil.h"
#include "external/TrackCovariance/SolTrack.h"

//
// Compare Kalman and full fit tracking errors

class KalmanCk{
	//
	// 
	// Author: F. Bedeschi, INFN-Pisa, Italy
	//  February 5, 2025
	//
private:
	//
	// Pointer to geometry
	SolGeom *fG;		// Geometry class
	Double_t fBz;		// Magnetic field		
	//
	// Resolution and multiple scattering activation flags
	Bool_t fRes;		// Detector resolution flag
	Bool_t fMS;		// Multiple scattering flag
	//
	// Option to use old covariance calculation
	Bool_t fOld;	// If true select old method
public:
	//
	// Constructors
	KalmanCk(SolGeom *G);
	// Destructor
	~KalmanCk();
	//
	// arrays for pt scans
	Int_t fNpt_Fa;			// Nr. of pt points in pt scans
	std::vector<Double_t> fPt_fixA;	// Array of pt points in pt scans
	Int_t fNangFa;			// Nr. of angles in pt scans
	std::vector<Double_t> fAngFa;	// Array of angles in pt scans
	// Arrays for angle scans
	Int_t fNang_Fpt;		// Nr. of angle points in angle scans
	std::vector<Double_t>fAng_fixPt;// Array of angles in angle scans
	Int_t fNptFpt;			// Nr. of pt angle scans
	std::vector<Double_t> fPtFpt;	// Array of pt in angle scans
	//
	// Standard resolution graphs for pt scans
	std::vector<TGraph*> gs_D_Pt;	// D resolution graphs for pt scans
	std::vector<TGraph*> gs_Phi0_Pt;// Phi0 resolution graphs for pt scans
	std::vector<TGraph*> gs_Pt_Pt;	// Pt resolution graphs for pt scans
	std::vector<TGraph*> gs_z0_Pt;	// z0 resolution graphs for pt scans
	std::vector<TGraph*> gs_cot_Pt;	// cot(theta) graphs for pt scans
	//
	// Kalman resolution graphs for pt scans
	std::vector<TGraph*> gk_D_Pt;	// D resolution graphs for pt scans
	std::vector<TGraph*> gk_Phi0_Pt;// Phi0 resolution graphs for pt scans
	std::vector<TGraph*> gk_Pt_Pt;	// Pt resolution graphs for pt scans
	std::vector<TGraph*> gk_z0_Pt;	// z0 resolution graphs for pt scans
	std::vector<TGraph*> gk_cot_Pt;	// cot(theta) graphs for pt scans
	//
	// Standard resolution graphs for angle scans
	std::vector<TGraph*> gs_D_Ang;		// D resolution graphs for angle scans
	std::vector<TGraph*> gs_Phi0_Ang;	// Phi0 resolution graphs for angle scans
	std::vector<TGraph*> gs_Pt_Ang;		// Pt resolution graphs for angle scans
	std::vector<TGraph*> gs_z0_Ang;		// z0 resolution graphs for angle scans
	std::vector<TGraph*> gs_cot_Ang;	// cot(theta) graphs for angle scans
	//
	// Kalman resolution graphs for angle scans
	std::vector<TGraph*> gk_D_Ang;		// D resolution graphs for angle scans
	std::vector<TGraph*> gk_Phi0_Ang;	// Phi0 resolution graphs for angle scans
	std::vector<TGraph*> gk_Pt_Ang;		// Pt resolution graphs for angle scans
	std::vector<TGraph*> gk_z0_Ang;		// z0 resolution graphs for angle scans
	std::vector<TGraph*> gk_cot_Ang;	// cot(theta) graphs for angle scans
	//
	//
	// Methods
	//
	//
	// Utility
	void GrSetup(TGraph* &g, Double_t yMax, Int_t Color, TString title, TString axis);
	//
	void SetMode(Bool_t Res, Bool_t MS);	// Set detector resolution and/or multiple scatt.
	void SetOld(Bool_t Old){ fOld = Old;};	// Select old Covariance calculation
											// Call before Fill !!! Default is true, true
	void DrawPtScan(Int_t i);	// Draw tracks with fixed pt = fPt_fixA[i];
	void Fill();				// Fill resolution graphs
	void Print();				// Display all graphs
};

#endif
