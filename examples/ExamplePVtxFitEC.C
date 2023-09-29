/*
Example of using vertex fitter class to fit primary vertex
assumed to be generated with Pythia8/ee_zh_smr-shf.cmn
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)

#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "modules/TrackCovariance.h"
#include "external/TrackCovariance/TrkUtil.h"
#include "external/TrackCovariance/VertexFit.h"
#include "examplles/classes/VtxPulls.h"

#endif



//------------------------------------------------------------------------------

void ExamplePVtxFitEC(const char* inputFile, Int_t Nevent = 5)
{
	//
	// Beam constraint 
	// (consistent with generation with ee_zh_smr-shf.cmnd)
	//
	// Mean beam position
	TVectorD xpvc(3);
	xpvc(0) = 1.0;
	xpvc(1) = -2.0;
	xpvc(2) = 10.0;
	// Interaction region covariance
	TMatrixDSym covpvc(3); covpvc.Zero();
	covpvc(0, 0) = 0.0097 * 0.0097;
	covpvc(1, 1) = 2.55e-05 * 2.55e-05;
	covpvc(2, 2) = 0.64 * 0.64;

	//
	// Create chain of root trees
	TChain chain("Delphes");
	chain.Add(inputFile);

	// Create object of class ExRootTreeReader
	ExRootTreeReader* treeReader = new ExRootTreeReader(&chain);
	Long64_t numberOfEntries = treeReader->GetEntries();

	// Get pointers to branches used in this analysis
	TClonesArray* branchGenPart = treeReader->UseBranch("Particle");
	TClonesArray* branchTrack = treeReader->UseBranch("Track");

	// Book histograms
	VtxPulls* Pulls = new VtxPulls();
	//
	// Loop over all events
	Int_t Nev = TMath::Min(Nevent, (Int_t)numberOfEntries);
	for (Int_t entry = 0; entry < Nev; ++entry)
	{
		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);
		// 
		Int_t Ntr = 0;	// # of tracks from primary vertex
		Int_t NtrG = branchTrack->GetEntries();
		TVectorD** pr = new TVectorD * [NtrG];
		TMatrixDSym** cv = new TMatrixDSym * [NtrG];
		TVectorD vTrue(3);
		std::vector<TVectorD> pGen;	// True momenta
		//
		// Print event numbers
		if(entry % 500 == 0){
			std::cout<<std::endl;
			std::cout <<"Event " << entry << std::endl;
			std::cout<<"Nr. of tracks " << NtrG << std::endl;
		}
		//
		// True vertex
		Double_t xpv, ypv, zpv;
		// If event contains at least 1 track
		//
		if (branchTrack->GetEntries() > 0)
		{
			// Loop on tracks
			for (Int_t it = 0; it < branchTrack->GetEntries(); it++)
			{
				Track* trk = (Track*)branchTrack->At(it);
				//
				// Get associated generated particle
				GenParticle* gp = (GenParticle*)trk->Particle.GetObject();
				TVector3 xg(1.e-3*gp->X,1.e-3*gp->Y,1.e-3*gp->Z); 	// mm -> meters
				TVector3 pg(gp->Px,gp->Py,gp->Pz);
				Double_t Q = (Double_t)gp->Charge;
				Double_t Bz = 2.0;
				TVectorD genParM =TrkUtil:: XPtoPar(xg, pg, Q, Bz);
				TVectorD genPar = TrkUtil::ParToMm(genParM); 		// -> back to mm

				// 
				// Check if original track is from primary vertex
				//
				// Position of origin in mm
				Double_t x = gp->X;
				Double_t y = gp->Y;
				Double_t z = gp->Z;
				Bool_t prim = kTRUE;	// Is primary?
				Int_t mp = gp->M1;	// Mother
				while(mp>0){
					GenParticle* gm = (GenParticle*)branchGenPart->At(mp);
					Double_t xm = gm->X;
                                	Double_t ym = gm->Y;
                                	Double_t zm = gm->Z;
					if(x!=xm || y!=ym || z!=zm){
						prim = kFALSE;	// Not primary
						break;
					}else mp = gm->M1;
				}

				//
				// group tracks originating from the primary vertex
				if (prim)
				{
					//
					// Store true primary vertex for this event
					xpv = x;
					ypv = y;
					zpv = z; 
					vTrue(0) = x;
					vTrue(1) = y;
					vTrue(2) = z;
					//
					// Reconstructed track parameters
					Double_t obsD0 = trk->D0;
					Double_t obsPhi = trk->Phi;
					Double_t obsC = trk->C;
					Double_t obsZ0 = trk->DZ;
					Double_t obsCtg = trk->CtgTheta;
					Double_t oPar[5] = { obsD0, obsPhi, obsC, obsZ0, obsCtg };
					TVectorD obsPar(5, oPar);	// Fill observed parameters
					//
					pr[Ntr] = new TVectorD(obsPar);
					cv[Ntr] = new TMatrixDSym(trk->CovarianceMatrix());
					Ntr++;
					// Store true momenta
					GenParticle* gen = (GenParticle*) trk->Particle.GetObject();
					TVectorD GenP(3);
					GenP(0) = gen->Px;
					GenP(1) = gen->Py;
					GenP(2) = gen->Pz;
					pGen.push_back(GenP);
				}
			}		// End loop on tracks
		}
		if(entry % 100 == 0)std::cout<<"Event # "<<entry<<std::endl;
		//
		// Fit primary vertex with beam constraint
		//
		Int_t MinTrk = 1;	// Minumum # tracks for vertex fit
		if (Ntr >= MinTrk) {
			VertexFit* Vtx = new VertexFit(Ntr, pr, cv);
			Vtx->AddVtxConstraint(xpvc, covpvc);
			TVectorD xvtx = Vtx->GetVtx();
			//std::cout<<"N tracks = "<<Ntr<<", X,Y,Z in: "<<xpvc(0)<<", "<<xpvc(1)<<", "<<xpvc(2)<<std::endl;
			//std::cout<<"X, Y, Z out: "<<xvtx(0)<<", "<<xvtx(1)<<", "<<xvtx(2)<<std::endl;
			for(Int_t i=0; i<Ntr; i++){
				TVectorD pMom = pGen[i];
				//pMom.Print();
			}
			Pulls->Fill(vTrue, pGen, Vtx);
			delete Vtx;
		}

		//
		// Cleanup
		for (Int_t i = 0; i < Ntr; i++) delete pr[i];
		for (Int_t i = 0; i < Ntr; i++) delete cv[i];
		delete[] pr;
		delete[] cv;
		pGen.clear();
	}

	//
	// Show resulting histograms
	//
	Pulls->Print();
}
