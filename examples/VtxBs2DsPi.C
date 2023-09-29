/*
Simple macro showing how to access branches from the delphes output root file,
loop over events, and plot track pulls and geometrical acceptance.

root -l examples/Example6.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include <TClonesArray.h>
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "modules/TrackCovariance.h"
#include "external/TrackCovariance/TrkUtil.h"
#include "external/TrackCovariance/VertexFit.h"
#include "external/TrackCovariance/VertexMore.h"
#include "examples/classes/VState.h"
#include "examples/classes/VList.h"
#include "examples/classes/BsPulls.h"
#endif
//
// Global variables
//
// Particle codes
	Int_t PIpID = 211; 	Int_t PImID = -211;
	Int_t KpID  = 321; 	Int_t KmID  = -321;
	Int_t PhiID = 333;
	Int_t DSpID = 431;	Int_t DSmID = -431;
	Int_t BS0ID = 531; 	Int_t BSbID = -531;

//
// Find Bs mesons decaying ro Ds and charged pion
std::vector<VList*> FindBsToDsPi(TClonesArray* branchTrack, TClonesArray* branchGenPart, Int_t prtOpt)
{
	//
	// Fill Phi list
	std::vector<Int_t> tPhiPID;
	tPhiPID.push_back(KpID); tPhiPID.push_back(KmID);
	std::vector<VList*> vPhiLink;
	VList* vPhi = new VList(PhiID, tPhiPID, vPhiLink, branchTrack, branchGenPart);
	//
	// Fill Ds+ list
	std::vector<Int_t> tDspPID;
	tDspPID.push_back(PIpID);
	std::vector<VList*> vDspLink;
	vDspLink.push_back(vPhi);
	VList* vDsp = new VList(DSpID, tDspPID, vDspLink, branchTrack, branchGenPart);
	//
	// Fill Ds- list
	std::vector<Int_t> tDsmPID;
	tDsmPID.push_back(PImID);
	std::vector<VList*> vDsmLink;
	vDsmLink.push_back(vPhi);
	VList* vDsm = new VList(DSmID, tDsmPID, vDsmLink, branchTrack, branchGenPart);
	//
	// All Bs lists
	std::vector<VList*>vBsAll;
	//
	// Fill Bs list
	std::vector<Int_t> tBsPID;
	tBsPID.push_back(PIpID);
	std::vector<VList*> vBsLink;
	vBsLink.push_back(vDsm);
	VList* vBs0 = new VList(BS0ID, tBsPID, vBsLink, branchTrack, branchGenPart);
	if(vBs0->GetNvtx() > 0)vBsAll.push_back(vBs0);	
	//
	// Fill Bs mixed list
	std::vector<Int_t> tBsPIDm;
	tBsPIDm.push_back(PImID);
	std::vector<VList*> vBsLinkm;
	vBsLinkm.push_back(vDsp);
	VList* vBs0m = new VList(BS0ID, tBsPIDm, vBsLinkm, branchTrack, branchGenPart);
	if(vBs0m->GetNvtx() > 0)vBsAll.push_back(vBs0m);
	//
	// Fill Bs-bar list
	std::vector<Int_t> tBsbPID;
	tBsbPID.push_back(PImID);
	std::vector<VList*> vBsbLink;
	vBsbLink.push_back(vDsp);
	VList* vBs0b = new VList(BSbID, tBsbPID, vBsbLink, branchTrack, branchGenPart);
	if(vBs0b->GetNvtx() > 0)vBsAll.push_back(vBs0b);
	//
	// Fill Bs-bar mixed list
	std::vector<Int_t> tBsbPIDm;
	tBsbPIDm.push_back(PIpID);
	std::vector<VList*> vBsbLinkm;
	vBsbLinkm.push_back(vDsm);
	VList* vBs0bm = new VList(BSbID, tBsbPIDm, vBsbLinkm, branchTrack, branchGenPart);
	if(vBs0bm->GetNvtx() > 0) vBsAll.push_back(vBs0bm);
	
	//
	// Printout
	//
	if(prtOpt == 2){
		// Phi
		Int_t nPhi = vPhi->GetNvtx();
		if(nPhi>0){
			std::cout << "Number of Phi found: " << nPhi << std::endl;
			for (Int_t i = 0; i < nPhi; i++)vPhi->GetVState(i)->Print();
		}
		// Ds+
		Int_t nDsp = vDsp->GetNvtx();
		if(nDsp>0){
			std::cout << "Number of Ds+ found: " << nDsp << std::endl;
			for (Int_t i = 0; i < nDsp; i++)vDsp->GetVState(i)->Print();
		}
		// Ds-
		Int_t nDsm = vDsm->GetNvtx();
		if(nDsm>0){
			std::cout << "Number of Ds- found: " << nDsm << std::endl;
			for (Int_t i = 0; i < nDsm; i++)vDsm->GetVState(i)->Print();
		}
	}
	if(prtOpt >=1){
		// Bs0
		Int_t nBs0 = vBs0->GetNvtx();
		if(nBs0>0){
			std::cout << "Number of Bs0 found: " << nBs0 << std::endl;
			for (Int_t i = 0; i < nBs0; i++)vBs0->GetVState(i)->Print();
		}
		// Bs0 mixed
		Int_t nBs0m = vBs0m->GetNvtx();
		if(nBs0m>0){
			std::cout << "Number of mixed Bs0 found: " << nBs0m << std::endl;
			for (Int_t i = 0; i < nBs0m; i++)vBs0m->GetVState(i)->Print();
		}
		// Bs0-bar
		Int_t nBs0b = vBs0b->GetNvtx();
		if(nBs0b>0){
			std::cout << "Number of Bs0b found: " << nBs0b << std::endl;
			for (Int_t i = 0; i < nBs0b; i++)vBs0b->GetVState(i)->Print();
		}
		// Bs0-bar mixed
		Int_t nBs0bm = vBs0bm->GetNvtx();
		if(nBs0bm>0){
			std::cout << "Number of mixed Bs0b found: " << nBs0bm << std::endl;
			for (Int_t i = 0; i < nBs0bm; i++)vBs0bm->GetVState(i)->Print();
		}
	}
//
	return vBsAll;
}

//
// Track to vector & cov
void TrkToVector(Track* tr, TVectorD &par, TMatrixDSym &cov)
{
	Double_t trPar[5] = {tr->D0,tr->Phi,tr->C,tr->DZ,tr->CtgTheta};
	TVectorD tPar(5,trPar); 
	par = tPar;
	cov = tr->CovarianceMatrix();
}
//------------------------------------------------------------------------------

void VtxBs2DsPi(const char* inputFile, Int_t Nevent = 10, Int_t prtOpt = 0)
{
	gSystem->Load("libDelphes");

	// Create chain of root trees
	TChain chain("Delphes");
	chain.Add(inputFile);

	// Create object of class ExRootTreeReader
	ExRootTreeReader* treeReader = new ExRootTreeReader(&chain);
	Long64_t numberOfEntries = treeReader->GetEntries();

	// Get pointers to branches used in this analysis
	TClonesArray* branchGenPart = treeReader->UseBranch("Particle");
	TClonesArray* branchTrack = treeReader->UseBranch("Track");

	// 
	// Initialize histograms
	BsPulls* BsHist = new BsPulls();

	//
	// Loop over all events
	Int_t Nev = TMath::Min(Nevent, (Int_t)numberOfEntries);
	for (Int_t entry = 0; entry < Nev; ++entry)
	{
		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);
		Int_t NtrG = branchTrack->GetEntries();
		std::cout<<std::endl;
		std::cout <<"Event " << entry << std::endl;
		//std::cout<<"Nr. of tracks " << NtrG << std::endl;
		//
		// test Particle branch
		Int_t Ngen = branchGenPart->GetEntries();
		//std::cout << "Nr. of generated particles: " << Ngen << std::endl;

		// 
		// Find Bs-> Ds pi including mixing
		std::vector<VList*>vBsAll = FindBsToDsPi(branchTrack, branchGenPart, prtOpt);

		//
		// Chain fit Bs --> Ds pi
		for(UInt_t i=0; i<vBsAll.size(); i++){	// Loop on Bs types
			VList* vBs = vBsAll[i];		// List of found Bs of this type
			Int_t nBs = vBs->GetNvtx();	// Nr. found
			for(Int_t j=0; j<nBs; j++){	// Loop on found Bs in event
				//
				// Extract tracks for fitting
				const Int_t nDsT = 3;	// Tracks in Ds fit
				Track** tDs = new Track*[nDsT];
				GenParticle** pDs = new GenParticle*[nDsT];
				const Int_t nBsT = 2;	// Tracks in Bs fit
				Track** tBs = new Track*[nBsT];
				GenParticle** pBs = new GenParticle*[nBsT];
				//
				VState* DsPiState   = vBs->GetVState(j);	// Specific final state Ds pi
				Track*  tBsPi       = DsPiState ->GetTrk(0);	// Pion from Bs
				tBs[0] = tBsPi;
				GenParticle* pBsPi  = DsPiState ->GetGen(0);	// Associated generated particle
				pBs[0] = pBsPi;
				VState* PhiPiState  = DsPiState ->GetVState(0);	// Ds from Bs
				GenParticle* pBsDs  = PhiPiState->GetMother();
				pBs[1] = pBsDs;
				Track*  tDsPi       = PhiPiState->GetTrk(0);	// Pion from Ds
				tDs[0] = tDsPi;
				GenParticle* pDsPi  = PhiPiState->GetGen(0);	// Associated generated particle
				pDs[0] = pDsPi;
				VState* KpKmState   = PhiPiState->GetVState(0);	// Phi from Ds
				Track*  tPhiK1	    = KpKmState ->GetTrk(0);	// First  K from Phi
				tDs[1] = tPhiK1;
				GenParticle* pPhiK1 = KpKmState ->GetGen(0);	// Associated generated particle
				pDs[1] = pPhiK1;
				Track*  tPhiK2	    = KpKmState ->GetTrk(1);	// Second K from Phi
				tDs[2] = tPhiK2;
				GenParticle* pPhiK2 = KpKmState ->GetGen(1);	// Associated generated particle
				pDs[2] = pPhiK2;
				//
				// Load Ds tracks
				TVectorD* tDsPar[nDsT];
				TMatrixDSym* tDsCov[nDsT];
				for(Int_t k=0; k<nDsT; k++){
					TVectorD par(5); TMatrixDSym cov(5);
					TrkToVector(tDs[k], par, cov);
					tDsPar[k] = new TVectorD(par);
					tDsCov[k] = new TMatrixDSym(cov);
				}
				// Fit Ds vertex
				VertexFit* vDs = new VertexFit(nDsT, tDsPar, tDsCov);
				Double_t DsChi2 = vDs->GetVtxChi2();		// Ds fit Chi2
				// More fitting
				Bool_t Units = kTRUE;		// Set to mm		
				VertexMore* VMDs = new VertexMore(vDs,Units);
				// Mass constraint if requested
				Bool_t MCst = kTRUE;			// Mass constraint flag
				if(MCst){
					Double_t DsMass = pBsDs->Mass;	// Ds mass
					Double_t DsMasses[nDsT]; Int_t DsList[nDsT];
					for(Int_t k=0; k<nDsT; k++){
						DsMasses[k] = pDs[k]->Mass;
						DsList[k]   = k;
					}
					VMDs->AddMassConstraint(DsMass,  nDsT,  DsMasses,  DsList);
					VMDs->MassConstrFit();
				}
				TVectorD rDv = VMDs->GetXv();			// Ds vertex
				Double_t RDs = TMath::Sqrt(rDv[0]*rDv[0]+rDv[1]*rDv[1]);

				//
				// Load Bs tracks
				TVectorD* tBsPar[nBsT];
				TMatrixDSym* tBsCov[nBsT];
				TVectorD par(5); TMatrixDSym cov(5);
				TrkToVector(tBs[0], par, cov);
				tBsPar[0] = new TVectorD(par);			// Bs pion	
				tBsCov[0] = new TMatrixDSym(cov);
				tBsPar[1] = new TVectorD(VMDs->GetVpar());	// Ds from previous fit
				tBsCov[1] = new TMatrixDSym(VMDs->GetVcov());
				//
				// Fit Bs vertex
				VertexFit* vBs = new VertexFit(nBsT, tBsPar, tBsCov);
;				vBs->SetStartR(RDs);
				Double_t BsChi2 = vBs->GetVtxChi2();
				// More fitting
				VertexMore* VMBs = new VertexMore(vBs,Units);

				//
				// Fill histograms
				BsHist->Fill(DsPiState, VMDs, VMBs);
				
			} 	// End loop on found Bs	
		} 	// End loop on Bs types

	}	// End event loop

//
// Print histograms
	BsHist->Print();
}
