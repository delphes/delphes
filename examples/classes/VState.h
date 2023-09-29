//
#ifndef G__VSTATE_H
#define G__VSTATE_H
//
#include <vector>
#include <iostream>
#include <TString.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include "classes/DelphesClasses.h"
//
// Class to store tracks belonging to a common vertex

class VState{
	//
	// 
	// Author: F. Bedeschi, INFN-Pisa, Italy
	//  October 25, 2022
	//
private:
	Int_t fMotherPID;
	GenParticle* fMin;
	std::vector<Track*> ftLink;
	std::vector<GenParticle*> fpLink;
	std::vector<VState*> fvLink;
	//	
public:
	//
	// Constructors
	VState(GenParticle* gMin, std::vector<Track*>& tLink, std::vector<VState*>& vLink);
	// Destructor
	~VState();
	//
	//
	// Methods
	GenParticle* GetMother() { return fMin; };		// Return mother pointer
	Int_t GetMotherPID() { return fMotherPID; };		// Return mother PID
	Int_t GetNtrk() { return ftLink.size(); };		// Return # tracks in vertex
	Track* GetTrk(Int_t i) { return ftLink[i]; };		// Return pointer to track
	GenParticle* GetGen(Int_t i) { return fpLink[i]; };	// Return pointer to generated particle
	Int_t GetNvtx() { return fvLink.size(); };		// Return # vertices in vertex
	VState* GetVState(Int_t i) { return fvLink[i]; };	// Return pointer to vertex
	//
	static TString GetPDGname(Int_t PID);			// Get Particle name
	void Print();						// Printout this state
};

#endif
