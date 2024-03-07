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
#include "examples/classes/B0KsKsPulls.h"
#include "examples/classes/KsPulls.h"
#endif
//
// Global variables
//
// Particle codes
	Int_t PIpID = 211; 	Int_t PImID = -211;
	Int_t KsID  = 310; 
	Int_t B0ID = 511; 	Int_t B0bID = -511;

//
// Find Bs mesons decaying ro Ds and charged pion
std::vector<VList*> FindB0ToKsKs(TClonesArray* branchTrack, TClonesArray* branchGenPart, Int_t prtOpt)
{
	//
	// Fill Ks list
	std::vector<Int_t> tKsPID;
	tKsPID.push_back(PIpID); tKsPID.push_back(PImID);
	std::vector<VList*> vKsLink;
	VList* vKs = new VList(KsID, tKsPID, vKsLink, branchTrack, branchGenPart);
	//
	// All B0 lists
	std::vector<VList*>vB0All;
	//
	// Fill B0 list
	std::vector<Int_t> tB0PID;
	std::vector<VList*> vB0Link;
	vB0Link.push_back(vKs);
	vB0Link.push_back(vKs);
	VList* vB0 = new VList(B0ID, tB0PID, vB0Link, branchTrack, branchGenPart);
	if(vB0->GetNvtx() > 0)vB0All.push_back(vB0);
	//
	// Fill Bs-bar list
	std::vector<Int_t> tB0bPID;
	std::vector<VList*> vB0bLink;
	vB0bLink.push_back(vKs);
	vB0bLink.push_back(vKs);
	VList* vB0b = new VList(B0bID, tB0bPID, vB0bLink, branchTrack, branchGenPart);
	if(vB0b->GetNvtx() > 0)vB0All.push_back(vB0b);
	//
	// Printout
	//
	if(prtOpt == 2){
		// Ks
		Int_t nKs = vKs->GetNvtx();
		if(nKs>0){
			std::cout << "Number of Ks found: " << nKs << std::endl;
			for (Int_t i = 0; i < nKs; i++)vKs->GetVState(i)->Print();
		}
	}
	if(prtOpt >=1){
		// B0
		Int_t nB0 = vB0->GetNvtx();
		if(nB0>0){
			std::cout << "Number of B0 found: " << nB0 << std::endl;
			for (Int_t i = 0; i < nB0; i++)vB0->GetVState(i)->Print();
		}
		// B0-bar
		Int_t nB0b = vB0b->GetNvtx();
		if(nB0b>0){
			std::cout << "Number of B0b found: " << nB0b << std::endl;
			for (Int_t i = 0; i < nB0b; i++)vB0b->GetVState(i)->Print();
		}
	}
//
	return vB0All;
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

void VtxB02KsKs(const char* inputFile, Int_t Nevent = 10, Int_t prtOpt = 0)
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
	B0KsKsPulls* B0Hist = new B0KsKsPulls();
	KsPulls* KsPlots = new KsPulls();
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
		// Particle branch
		Int_t Ngen = branchGenPart->GetEntries();
		//std::cout << "Nr. of generated particles: " << Ngen << std::endl;

		// 
		// Find B-> Ks Ks 
		std::vector<VList*>vB0All = FindB0ToKsKs(branchTrack, branchGenPart, prtOpt);
		//
		// Chain fit B0 --> Ks Ks
		for(UInt_t i=0; i<vB0All.size(); i++){	// Loop on found B0 types
			VList* vB0List = vB0All[i];	// List of found B0 of this type
			Int_t nB0 = vB0List->GetNvtx();	// Nr. found
			for(Int_t j=0; j<nB0; j++){	// Loop on found B0 in event
				//
				// Extract tracks for fitting
				const Int_t nKsT = 2;			// Tracks in Ks fit
				Track* tKs1[nKsT];			// First Ks
				GenParticle* pKs1[nKsT];
				Track* tKs2[nKsT];			// Second Ks
				GenParticle* pKs2[nKsT];
				const Int_t nB0T = 2;			// Tracks in B0 fit (neutral)
				GenParticle* pB0[nB0T];
				//
				VState* B0State     = vB0List->GetVState(j);	// Specific final state KsKs
				// 1st Ks
				VState* Ks1State    = B0State->GetVState(0);	// 1st Ks
				pB0[0]              = Ks1State->GetMother();	// Pointer to Ks genparticle
				tKs1[0] 	    = Ks1State->GetTrk(0);	// 1st pion associated track
				pKs1[0] 	    = Ks1State->GetGen(0);	// Associated generated particle
				tKs1[1]	   	    = Ks1State->GetTrk(1);	// 2nd pion associated track
				pKs1[1] 	    = Ks1State->GetGen(1);	// Associated generated particle	
				// 2nd Ks
				VState* Ks2State    = B0State->GetVState(1);	// 2nd Ks
				pB0[1]              = Ks2State->GetMother();	// Pointer to Ks genparticle
				tKs2[0] 	    = Ks2State->GetTrk(0);	// 1st pion associated track
				pKs2[0] 	    = Ks2State->GetGen(0);	// Associated generated particle
				tKs2[1]	   	    = Ks2State->GetTrk(1);	// 2nd pion associated track
				pKs2[1] 	    = Ks2State->GetGen(1);	// Associated generated particle
				//
				// Load 1st Ks tracks
				TVectorD* tKs1Par[nKsT];
				TMatrixDSym* tKs1Cov[nKsT];
				for(Int_t k=0; k<nKsT; k++){
					TVectorD par(5); TMatrixDSym cov(5);
					TrkToVector(tKs1[k], par, cov);
					tKs1Par[k] = new TVectorD(par);
					tKs1Cov[k] = new TMatrixDSym(cov);
				}
				// Fit 1st Ks vertex
				VertexFit* vKs1 = new VertexFit(nKsT, tKs1Par, tKs1Cov);
				Double_t R1 = TMath::Sqrt(pow(tKs1[0]->XFirstHit,2)+pow(tKs1[0]->YFirstHit,2));
				Double_t R2 = TMath::Sqrt(pow(tKs1[1]->XFirstHit,2)+pow(tKs1[1]->YFirstHit,2));
				Double_t Rmin = TMath::Min(R1,R2);		// Smallest radius of first hit
				Rmin = TMath::Sqrt(pow(pKs1[0]->X,2)+pow(pKs1[0]->Y,2));	// Use real vertex
				vKs1->SetStartR(Rmin);				// Protect against double crossings
				Double_t Ks1Chi2 = vKs1->GetVtxChi2();		// Ks fit Chi2
				//std::cout<<"Ks 1 chi2 = "<<Ks1Chi2<<std::endl;
				// More fitting
				Bool_t Units = kTRUE;		// Set to mm		
				VertexMore* VMKs1 = new VertexMore(vKs1,Units);
				// Mass constraint if requested
				Bool_t MCst = kTRUE;			// Mass constraint flag
				if(MCst){
					Double_t KsMass = pB0[0]->Mass;	// Ks mass
					Double_t KsMasses[nKsT]; Int_t KsList[nKsT];
					for(Int_t k=0; k<nKsT; k++){
						KsMasses[k] = pKs1[k]->Mass;
						KsList[k]   = k;
					}
					//std::cout<<"K0S nr 1 mass= "<<KsMass
					//<<" Pi1= "<<KsMasses[0]<<" Pi2= "<<KsMasses[1]<<std::endl;
					VMKs1->AddMassConstraint(KsMass,  nKsT,  KsMasses,  KsList);
					VMKs1->MassConstrFit();
				}
				TVectorD rKv1 = VMKs1->GetXv();			// Ks vertex
				Double_t RKs1 = TMath::Sqrt(rKv1[0]*rKv1[0]+rKv1[1]*rKv1[1]);
				KsPlots->Fill(Ks1State, VMKs1);
				//std::cout<<"First Ks fit completed"<<std::endl;
				//
				// Load 2nd Ks tracks
				TVectorD* tKs2Par[nKsT];
				TMatrixDSym* tKs2Cov[nKsT];
				for(Int_t k=0; k<nKsT; k++){
					TVectorD par(5); TMatrixDSym cov(5);
					TrkToVector(tKs2[k], par, cov);
					tKs2Par[k] = new TVectorD(par);
					tKs2Cov[k] = new TMatrixDSym(cov);
				}
				// Fit 2nd Ks vertex
				VertexFit* vKs2 = new VertexFit(nKsT, tKs2Par, tKs2Cov);
				R1 = TMath::Sqrt(pow(tKs2[0]->XFirstHit,2)+pow(tKs2[0]->YFirstHit,2));
				R2 = TMath::Sqrt(pow(tKs2[1]->XFirstHit,2)+pow(tKs2[1]->YFirstHit,2));
				Rmin = TMath::Min(R1,R2);		// Smallest radius of first hit
				Rmin = TMath::Sqrt(pow(pKs2[0]->X,2)+pow(pKs2[0]->Y,2));	// Use real vertex
				vKs2->SetStartR(Rmin);			// Protect against double crossings
				Double_t Ks2Chi2 = vKs2->GetVtxChi2();	// Ks fit Chi2
				//std::cout<<"Ks 2 chi2 = "<<Ks2Chi2<<std::endl;
				// More fitting
				Units = kTRUE;		// Set to mm		
				VertexMore* VMKs2 = new VertexMore(vKs2,Units);
				// Mass constraint if requested
				//
				if(MCst){
					Double_t KsMass = pB0[1]->Mass;	// Ks mass
					Double_t KsMasses[nKsT]; Int_t KsList[nKsT];
					for(Int_t k=0; k<nKsT; k++){
						KsMasses[k] = pKs2[k]->Mass;
						KsList[k]   = k;
					}
					//std::cout<<"K0S nr 2 mass= "<<KsMass
					//<<" Pi1= "<<KsMasses[0]<<" Pi2= "<<KsMasses[1]<<std::endl;
					VMKs2->AddMassConstraint(KsMass,  nKsT,  KsMasses,  KsList);
					VMKs2->MassConstrFit();
				}
				TVectorD rKv2 = VMKs2->GetXv();			// Ks vertex
				Double_t RKs2 = TMath::Sqrt(rKv2[0]*rKv2[0]+rKv2[1]*rKv2[1]);
				KsPlots->Fill(Ks2State, VMKs2);
				//std::cout<<"Second Ks fit completed"<<std::endl;
				//
				// Load B0 tracks
				TVectorD* tB0Par[nB0T];
				TMatrixDSym* tB0Cov[nB0T];
				tB0Par[0] = new TVectorD(VMKs1->GetVpar());	// 1st Ks from previous fit	
				tB0Cov[0] = new TMatrixDSym(VMKs1->GetVcov());
				tB0Par[1] = new TVectorD(VMKs2->GetVpar());	// 2nd Ks from previous fit
				tB0Cov[1] = new TMatrixDSym(VMKs2->GetVcov());
				//
				// Fit B0 vertex
				Bool_t Charged[nB0T] = {kFALSE, kFALSE};
				VertexFit* vB0 = new VertexFit(nB0T, tB0Par, tB0Cov, Charged);
				Double_t B0Chi2 = vB0->GetVtxChi2();
				TVectorD rvB0 = vB0->GetVtx();
				//std::cout<<"Rec: "<<rvB0(0)<<", "<<rvB0(1)<<", "<<rvB0(2)<<std::endl;
				// More fitting
				VertexMore* VMB0 = new VertexMore(vB0,Units);
				TVectorD rvmB0 = VMB0->GetXv();
				//std::cout<<"Rvm: "<<rvmB0(0)<<", "<<rvmB0(1)<<", "<<rvmB0(2)<<std::endl;

				//
				// Fill histograms
				B0Hist->Fill(B0State, VMKs1, VMKs2, VMB0);
			} 	// End loop on found B0	
		} 	// End loop on B0 types
	}	// End event loop

//
// Print histograms
	B0Hist->Print();
	KsPlots->Print();
}
