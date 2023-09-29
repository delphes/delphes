//
#ifndef G__VLIST_H
#define G__VLIST_H
//
#include <vector>
#include <iostream>
#include <TString.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include "VState.h"
//
// List all found vertices with given structure

class VList{
	//
	// 
	// Author: F. Bedeschi, INFN-Pisa, Italy
	//  October 25, 2022
	//
private:
	Int_t fMotherPID;			// Mother ID
	std::vector<VState*> fVtx;	// List of found final states with given structure
	Bool_t TrkMatch(GenParticle* gp_in, Track*& trk_out);		// Match generated particle to track
	TClonesArray* fTrack;
	TClonesArray* fGenPart;
	//	
public:
	//
	// Constructors
	VList(Int_t gmPID, std::vector<Int_t>& tPID, std::vector<VList*>& vLink, TClonesArray* branchTrack, TClonesArray* branchGenPart);
	// Destructor
	~VList();
	//
	//
	// Methods
	Int_t GetMotherPID() { return fMotherPID; };		// Return mother PID
	Int_t GetNvtx() { return fVtx.size(); };		// return number of vertices found
	VState* GetVState(Int_t i) { return fVtx[i]; };		// Return pointer to vertex
	//
};

#endif
