//
#ifndef G__KSPULLS_H
#define G__KSPULLS_H
//
#include <iostream>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include "external/TrackCovariance/VertexFit.h"
#include "external/TrackCovariance/VertexMore.h"
#include "examples/classes/VState.h"

//
// List all found vertices with given structure

class KsPulls{
	//
	// 
	// Author: F. Bedeschi, INFN-Pisa, Italy
	//  October 25, 2022
	//
private:
	// Pions from Ks after vertex fit
	// Parameters
	std::vector<TH1D*> h_Pi_Par;
	// Momentum
	TH1D* h_PiPx;
	TH1D* h_PiPy;
	TH1D* h_PiPz;
	// Correlations
	TH2D* h_PullPxCorrR;
	TH2D* h_PullPyCorrR;
	TH2D* h_PullPxCorrPt;
	TH2D* h_PullPyCorrPt;
	TH2D* h_PullPxCorrRD;
	TH2D* h_PullPyCorrRD;
	TH2D* h_PullPxCorrLm;
	TH2D* h_PullPyCorrLm;
	// Ks vertex distributions
	TH1D* h_KsRgen;
	TH1D* h_KsRrec;
	// Ks vertex pulls
	TH1D* h_KsXv;
	TH1D* h_KsYv;
	TH1D* h_KsZv;
	// Ks vertex momentum pulls
	TH1D* h_KsPx;
	TH1D* h_KsPy;
	TH1D* h_KsPz;
	// Ks parameters pulls
	std::vector<TH1D*> h_Ks_Par;
	// Ks mass
	TH1D* h_KsdGenMass;
	TH1D* h_KsdRecMass;
	TH1D* h_KsMassErr;
	TH1D* h_KsMassPull;
public:
	//
	// Constructors
	KsPulls();
	// Destructor
	~KsPulls();
	//
	//
	// Methods
	//
	void Fill(VState* KsState, VertexMore* vKs);
	void Print();
};

#endif
