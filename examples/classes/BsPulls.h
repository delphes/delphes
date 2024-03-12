//
#ifndef G__BSPULLS_H
#define G__BSPULLS_H
//
#include <iostream>
#include <TH1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include "external/TrackCovariance/TrkUtil.h"
#include "external/TrackCovariance/VertexFit.h"
#include "external/TrackCovariance/VertexMore.h"
#include "examples/classes/VState.h"

//
// List all found vertices with given structure

class BsPulls{
	//
	// 
	// Author: F. Bedeschi, INFN-Pisa, Italy
	//  October 25, 2022
	//
private:
	// Ds vertex
	TH1D* h_DsXv;
	TH1D* h_DsYv;
	TH1D* h_DsZv;
	// Ds vertex momentum
	TH1D* h_DsPx;
	TH1D* h_DsPy;
	TH1D* h_DsPz;
	// Ds tracks phases
	TH1D* h_DsTph;
	// Ds track parameters
	TH1D* h_DsParD;
	TH1D* h_DsParP0;
	TH1D* h_DsParC;
	TH1D* h_DsParZ0;
	TH1D* h_DsParCt;
	// Ds tracks momentum
	TH1D* h_DsPionPx;
	TH1D* h_DsPionPy;
	TH1D* h_DsPionPz;
	TH1D* h_DsKaonPx;
	TH1D* h_DsKaonPy;
	TH1D* h_DsKaonPz;
	// Bs vertex
	TH1D* h_BsXv;
	TH1D* h_BsYv;
	TH1D* h_BsZv;
	// Bs vertex momentum
	TH1D* h_BsPx;
	TH1D* h_BsPy;
	TH1D* h_BsPz;
	// Bs mass
	TH1D* h_BsMass;
	TH1D* h_BsMerr;
	TH1D* h_BsMpull;
public:
	//
	// Constructors
	BsPulls();
	// Destructor
	~BsPulls();
	//
	//
	// Methods
	//
	void Fill(VState* BsState, VertexMore* vDs, VertexMore* vBs);
	void Print();
};

#endif
