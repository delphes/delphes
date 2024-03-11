//
#ifndef G__B0KSKSPULLS_H
#define G__B0KSKSPULLS_H
//
#include <iostream>
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include "external/TrackCovariance/VertexFit.h"
#include "external/TrackCovariance/VertexMore.h"
#include "examples/classes/VState.h"

//
// List all found vertices with given structure

class B0KsKsPulls{
	//
	// 
	// Author: F. Bedeschi, INFN-Pisa, Italy
	//  October 25, 2022
	//
private:
	// 1st Ks vertex
	TH1D* h_Ks1Xv;
	TH1D* h_Ks1Yv;
	TH1D* h_Ks1Zv;
	// Ks vertex momentum
	TH1D* h_Ks1Px;
	TH1D* h_Ks1Py;
	TH1D* h_Ks1Pz;
	// 2nd Ks vertex
	TH1D* h_Ks2Xv;
	TH1D* h_Ks2Yv;
	TH1D* h_Ks2Zv;
	// Ks vertex momentum
	TH1D* h_Ks2Px;
	TH1D* h_Ks2Py;
	TH1D* h_Ks2Pz;
	// B0 vertex
	TH1D* h_B0Xv;
	TH1D* h_B0Yv;
	TH1D* h_B0Zv;
	// B0 vertex momentum
	TH1D* h_B0Px;
	TH1D* h_B0Py;
	TH1D* h_B0Pz;
	// B0 mass
	TH1D* h_B0Mass;
	TH1D* h_B0Merr;
	TH1D* h_B0Mpull;
	// B0 flight path
	TH1D* h_B0Lxy;
	TH1D* h_B0LxyS;
	TH1D* h_B0LxyPull;
	TH1D* h_B0Lxyz;
	TH1D* h_B0LxyzS;
	TH1D* h_B0LxyzPull;
	// Histogram list	
	TObjArray Hlist;
public:
	//
	// Constructors
	B0KsKsPulls();
	// Destructor
	~B0KsKsPulls();
	//
	//
	// Methods
	//
	void Fill(VState* BsState, VertexMore* vKs1, VertexMore* vKs2, VertexMore* vB0);
	void Print();
};

#endif
