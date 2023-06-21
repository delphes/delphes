//
#ifndef G__VTXPULLS_H
#define G__VTXPULLS_H
//
#include <iostream>
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include "external/TrackCovariance/VertexFit.h"
#include "external/TrackCovariance/VertexMore.h"

//
// Pulls of tracks parameters and vertex after vertex fit

class VtxPulls{
	//
	// 
	// Author: F. Bedeschi, INFN-Pisa, Italy
	//  March 6, 2023
	//
private:
	// x vertex
	TH1D* h_Xv;
	TH1D* h_Yv;
	TH1D* h_Zv;
	TH1D* h_Chi2;
	// Track parameters after fit
	TH1D* h_D;
	TH1D* h_P0;
	TH1D* h_C;
	TH1D* h_Z0;
	TH1D* h_Ct;
	// Track momenta after fit
	TH1D* h_Px;
	TH1D* h_Py;
	TH1D* h_Pz;
public:
	//
	// Constructors
	VtxPulls();
	// Destructor
	~VtxPulls();
	//
	//
	// Methods
	//
	void Fill(TVectorD vTrue, std::vector<TVectorD> pTrue, VertexFit *Vfit);
	void Print();
};

#endif
