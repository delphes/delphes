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
#include "examples/classes/KsPulls.h"
#endif
//
// Global variables
//
// Particle codes
	Int_t PIpID = 211; 	Int_t PImID = -211;
	Int_t KsID  = 310; 

//
// Find Ks-->pi pi
VList* FindKs(TClonesArray* branchTrack, TClonesArray* branchGenPart, Int_t prtOpt)
{
	//
	// Fill Ks list
	std::vector<Int_t> tKsPID;
	tKsPID.push_back(PIpID); tKsPID.push_back(PImID);
	std::vector<VList*> vKsLink;
	VList* vKs = new VList(KsID, tKsPID, vKsLink, branchTrack, branchGenPart);
	//
	//
	// Printout
	//
	if(prtOpt > 0){
		// Ks
		Int_t nKs = vKs->GetNvtx();
		if(nKs>0){
			std::cout << "Number of Ks found: " << nKs << std::endl;
			for (Int_t i = 0; i < nKs; i++)vKs->GetVState(i)->Print();
		}
	}
//
	return vKs;
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

void VtxKs(const char* inputFile, Int_t Nevent = 10, Int_t prtOpt = 0)
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
	KsPulls* KsHist = new KsPulls();

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
		// Find Ks 
		VList* vKsl = FindKs(branchTrack, branchGenPart, prtOpt);
		//
		// Fit all found Ks
		for(UInt_t i=0; i<vKsl->GetNvtx(); i++){	// Loop on found Ks types
			//
			// Extract tracks for fitting
			const Int_t nKsT = 2;			// Tracks in Ks decay
			Track* tKs[nKsT];			// Ks decay tracks
			GenParticle* pKs[nKsT];			// Generated Ks decay particles
			//
			VState* KsState = vKsl->GetVState(i);	// Specific final Ks state
			GenParticle* KsPart = KsState->GetMother();	// Pointer to Ks genparticle
			tKs[0] 		= KsState->GetTrk(0);	// 1st pion associated track
			pKs[0] 		= KsState->GetGen(0);	// Associated generated particle
			tKs[1]	   	= KsState->GetTrk(1);	// 2nd pion associated track
			pKs[1] 		= KsState->GetGen(1);	// Associated generated particle   
			//
			// Load Ks tracks
			TVectorD* tKsPar[nKsT];
			TMatrixDSym* tKsCov[nKsT];
			for(Int_t k=0; k<nKsT; k++){
				TVectorD par(5); TMatrixDSym cov(5);
				TrkToVector(tKs[k], par, cov);
				tKsPar[k] = new TVectorD(par);
				tKsCov[k] = new TMatrixDSym(cov);
			}
			// Fit Ks vertex
			VertexFit* vKs = new VertexFit(nKsT, tKsPar, tKsCov);
			Double_t R1 = TMath::Sqrt(pow(tKs[0]->XFirstHit,2)+pow(tKs[0]->YFirstHit,2));
			Double_t R2 = TMath::Sqrt(pow(tKs[1]->XFirstHit,2)+pow(tKs[1]->YFirstHit,2));
			Double_t Rmin = TMath::Min(R1,R2);		// Smallest radius of first hit
			Rmin = TMath::Sqrt(pow(pKs[0]->X,2)+pow(pKs[0]->Y,2));	// Use real vertex
			vKs->SetStartR(Rmin);				// Protect against double crossings
			Double_t KsChi2 = vKs->GetVtxChi2();		// Ks fit Chi2
			//std::cout<<"Ks  chi2 = "<<Ks1Chi2<<std::endl;
			// More fitting
			Bool_t Units = kTRUE;				// Set to mm		
			VertexMore* VMKs = new VertexMore(vKs,Units);
			// Mass constraint if requested
			Bool_t MCst = kTRUE;				// Mass constraint flag
			if(MCst){
				Double_t KsMass = KsPart->Mass;		// Ks mass
				Double_t KsMasses[nKsT]; Int_t KsList[nKsT];
				for(Int_t k=0; k<nKsT; k++){
					KsMasses[k] = pKs[k]->Mass;
					KsList[k]   = k;
				}
				//std::cout<<"K0S  mass= "<<KsMass
				//<<" Pi= "<<KsMasses[0]<<" Pi= "<<KsMasses[1]<<std::endl;
				VMKs->AddMassConstraint(KsMass,  nKsT,  KsMasses,  KsList);
				VMKs->MassConstrFit();
			}
				TVectorD rKv = VMKs->GetXv();			// Ks vertex
				Double_t RKs = TMath::Sqrt(rKv[0]*rKv[0]+rKv[1]*rKv[1]);
				//std::cout<<"Ks fit completed"<<std::endl;
				//
				//
				// Fill histograms
				KsHist->Fill(KsState, VMKs);
				
		} 	// End loop on found Ks	
	}	// End event loop

//
// Print histograms
	KsHist->Print();
}
