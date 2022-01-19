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
	Int_t Nbin = 100;
	// Vertex fit pulls
	TH1D* hXpull = new TH1D("hXpull", "Pull X vertex component", Nbin, -10., 10.);
	TH1D* hYpull = new TH1D("hYpull", "Pull Y vertex component", Nbin, -10., 10.);
	TH1D* hZpull = new TH1D("hZpull", "Pull Z vertex component", Nbin, -10., 10.);
	TH1D* hChi2 = new TH1D("hChi2", "Vertex #chi^{2}/N_{dof}", Nbin, 0., 10.);
	// Generation check
	TH1D* hXvpull = new TH1D("hXvpull", "Pull X generated vertex", Nbin, -10., 10.);
	TH1D* hYvpull = new TH1D("hYvpull", "Pull Y generated vertex", Nbin, -10., 10.);
	TH1D* hZvpull = new TH1D("hZvpull", "Pull Z generated vertex", Nbin, -10., 10.);

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
				}
			}		// End loop on tracks
		}
		//
		// Fit primary vertex with beam constraint
		//
		Int_t MinTrk = 3;	// Minumum # tracks for vertex fit
		if (Ntr >= MinTrk) {
			VertexFit* Vtx = new VertexFit(Ntr, pr, cv);
			Vtx->AddVtxConstraint(xpvc, covpvc);
			TVectorD xvtx = Vtx->GetVtx();
			TMatrixDSym covX = Vtx->GetVtxCov();
			Double_t Chi2 = Vtx->GetVtxChi2();
			Double_t Ndof = 2 * (Double_t)Ntr;
			delete Vtx;
			//
			Double_t PullX = (xvtx(0)-xpv) / TMath::Sqrt(covX(0, 0));
			Double_t PullY = (xvtx(1)-ypv) / TMath::Sqrt(covX(1, 1));
			Double_t PullZ = (xvtx(2)-zpv) / TMath::Sqrt(covX(2, 2));
			//
			Double_t PullXv = (xpvc(0)-xpv) / TMath::Sqrt(covpvc(0, 0));
			Double_t PullYv = (xpvc(1)-ypv) / TMath::Sqrt(covpvc(1, 1));
			Double_t PullZv = (xpvc(2)-zpv) / TMath::Sqrt(covpvc(2, 2));
			//
			//
			// Fill histograms
			hXpull->Fill(PullX);
			hYpull->Fill(PullY);
			hZpull->Fill(PullZ);
			hChi2->Fill(Chi2 / Ndof);
			//
			hXvpull->Fill(PullXv);
			hYvpull->Fill(PullYv);
			hZvpull->Fill(PullZv);
		}

		//
		// Cleanup
		for (Int_t i = 0; i < Ntr; i++) delete pr[i];
		for (Int_t i = 0; i < Ntr; i++) delete cv[i];
		delete[] pr;
		delete[] cv;
	}

	//
	// Show resulting histograms
	//
	TCanvas* Cnv = new TCanvas("Cnv", "Delphes primary vertex pulls", 50, 50, 900, 500);
	Cnv->Divide(2, 2);
	Cnv->cd(1); gPad->SetLogy(1); gStyle->SetOptFit(1111);
	hXpull->Fit("gaus"); hXpull->Draw();
	Cnv->cd(2); gPad->SetLogy(1); gStyle->SetOptFit(1111);
	hYpull->Fit("gaus"); hYpull->Draw();
	Cnv->cd(3); gPad->SetLogy(1); gStyle->SetOptFit(1111);
	hZpull->Fit("gaus"); hZpull->Draw();
	Cnv->cd(4); hChi2->Draw();
	//
	TCanvas* Cnv1 = new TCanvas("Cnv1", "Generated primary vertex pulls", 100, 100, 900, 500);
	Cnv1->Divide(3, 1);
	Cnv1->cd(1); gPad->SetLogy(1); gStyle->SetOptFit(1111);
	hXvpull->Fit("gaus"); hXvpull->Draw();
	Cnv1->cd(2); gPad->SetLogy(1); gStyle->SetOptFit(1111);
	hYvpull->Fit("gaus"); hYvpull->Draw();
	Cnv1->cd(3); gPad->SetLogy(1); gStyle->SetOptFit(1111);
	hZvpull->Fit("gaus"); hZvpull->Draw();
}
