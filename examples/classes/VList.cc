#include "VList.h"

//
// Vertex list configuration description
//
// Constructor
VList::VList(Int_t gmPID, std::vector<Int_t>& tPID, std::vector<VList*>& vLink, TClonesArray* branchTrack, TClonesArray* branchGenPart)
{
	fMotherPID = gmPID;		// Store mother pID
	fTrack = branchTrack;
	fGenPart = branchGenPart;
	//std::cout<<"VList: mother= "<< gmPID<<std::endl;
	//
	Int_t Ntot = tPID.size() + vLink.size();
	for (Int_t n = 0; n < branchGenPart->GetEntries(); n++) {
		GenParticle* gp = (GenParticle*)branchGenPart->At(n);
		if (gp->PID == gmPID) {						// Found mother
			Int_t Ndau = gp->D2 - gp->D1 + 1;		// # of daughters
			//std::cout<<"Found mother: Ndau= "<<Ndau<<", Ntot= "<<Ntot<<std::endl;
			std::vector<Bool_t>mStable(tPID.size());				// Matched stable particles flags
			std::vector<Track*>tStable(tPID.size());				// Matched stable tracks
			std::vector<Bool_t>mVertex(vLink.size());				// Matched vertices flag
			std::vector<VState*>tVertex(vLink.size());				// Matched vertices
			if (Ndau == Ntot) {
				//
				// Match daughters to input
				for (Int_t i = gp->D1; i <= gp->D2; i++) {
					GenParticle* gd = (GenParticle*)branchGenPart->At(i);
					// Match stable particles
					Track* trd = nullptr;
					if (tPID.size() >0 && gd->Status == 1){ 
						for (UInt_t j = 0; j < tPID.size(); j++) {
							if (gd->PID == tPID[j] && !mStable[j]) {
								//std::cout<<"Matched: "<<gd->PID<<std::endl;
								if (TrkMatch(gd, trd)) {
									mStable[j] = kTRUE;
									tStable[j] = trd;
									//std::cout<<"Matched+track: "<<gd->PID<<std::endl;
								}
							}
						}
					}
					// Match vertices
					if (vLink.size() > 0 && gd->Status != 1) {
						for (UInt_t j = 0; j < vLink.size(); j++) {
							if (gd->PID == vLink[j]->GetMotherPID() && !mVertex[j]) {
								Int_t Nvtx = vLink[j]->GetNvtx();
								//std::cout<<"Matched vtx: "<<gd->PID<<", N= "<<Nvtx<<std::endl;
								Bool_t Done = kFALSE;
								for(Int_t k=0; k<Nvtx; k++){
									VState* vs = vLink[j]->GetVState(k);
									//std::cout<<"vs.mother= "<<vs->GetMotherPID()<<std::endl;
									if(!Done && vs->GetMother() == gd){
										Bool_t NotUsed = kTRUE;
										for(UInt_t h=0; h<j; h++){
											if(tVertex[h]->GetMother() == gd)
											NotUsed = kFALSE;
										}
										if(NotUsed){
											mVertex[j] = kTRUE;
											tVertex[j] = vs;
											Done = kTRUE;
											//std::cout<<"Matched vtx/mother: "<<gd->PID<<std::endl;
										}
									}
								}
							}
						}

					}
				}
			}
			//
			// Check matching
			Bool_t tMatch = kTRUE;
			if (tPID.size() > 0) {
				for (UInt_t i = 0; i < tPID.size(); i++)if (!mStable[i])tMatch = kFALSE;
			}
			Bool_t vMatch = kTRUE;
			if (vLink.size() > 0) {
				for (UInt_t i = 0; i < vLink.size(); i++)if (!mVertex[i])vMatch = kFALSE;
			}
			Bool_t Matched = tMatch && vMatch;
			//
			// Fill new VState if matched
			if (Matched) {
				VState* vs = new VState(gp, tStable, tVertex);
				fVtx.push_back(vs);
			}
		}
	}
}

//
// Destructor
VList::~VList() {
	fVtx.clear();
}

//
// Match generated track to reconstructed track
//
Bool_t VList::TrkMatch(GenParticle* gp_in, Track* &trk_out)
{
	//
	// Find track (Trk) in track collection matching a given generated particle.
	// Returns kFALSE if not found.
	//
	Bool_t rCode = kFALSE;
	trk_out = nullptr;
	Int_t Ntr = fTrack->GetEntries();
	for (Int_t n = 0; n < Ntr; n++) {
		Track* trk = (Track*)fTrack->At(n);
		GenParticle* gp = (GenParticle*)trk->Particle.GetObject();
		if (gp == gp_in) {
			trk_out = trk;
			rCode = kTRUE;
			break;
		}
	}
	return rCode;
}
