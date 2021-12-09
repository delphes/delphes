//
#ifndef G__ACCEPTANCECLX_H
#define G__ACCEPTANCECLX_H
//
#include <TMath.h>
#include <TVectorF.h>
#include <TMatrixF.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TF2.h>
#include <iostream>
#include <vector>
#include "SolGeom.h"
#include "SolTrack.h"
//
// Class to create geometry for solenoid geometry

class AcceptanceClx {
	//
	// Class to handle storing and retieving of tracking acceptance
	//
private:
	TMatrixF fAcc;		// Acceptance matrix
	Int_t fNPtNodes;	// Numer of Pt nodes 
	TVectorF fPtArray;	// Array of Pt nodes
	Int_t fNThNodes;	// Numer of Theta nodes 
	TVectorF fThArray;	// Array of Theta nodes		(Theta in degrees)
	//
	// Service routines
	void VecInsert(Int_t i, Float_t x, TVectorF& Vec);
	void SplitPt(Int_t i, TVectorF &AccPt);	
	void SplitTh(Int_t i, TVectorF &AccTh);
public:
	//
	// Constructors
	AcceptanceClx(SolGeom *InGeo);				// Initialize arrays from geometry
	AcceptanceClx(TString InFile);				// Initialize from acceptance file
	// Destructor
	~AcceptanceClx();
	//
	// Accessors
	TMatrixF* GetAccMatrix() { return &fAcc; }
	Int_t GetNrPt() { return fNPtNodes; }
	Int_t GetNrTh() { return fNThNodes; }
	
	TVectorF* GetPtArray() { return &fPtArray; }
	TVectorF* GetThArray() { return &fThArray; }
	//
	// Read and write
	void ReadAcceptance(TString InFile);	// Stand alone usage
	void WriteAcceptance(TString OutFile);	// Stand alone usage
	void WriteAcceptance(TFile *OutFile);
	//
	// Function returning interpolated number of hit measurement layers
	Double_t HitNumber(Double_t pt, Double_t Theta);	// Theta in degrees
	Double_t HitNum(Double_t *x, Double_t *p);
};

#endif
