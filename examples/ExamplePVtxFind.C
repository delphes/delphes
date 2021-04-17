/*
Example of using vertex fitter class to fit primary vertex
assumed to be generated in (0,0,0)
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

void ExamplePVtxFind(const char* inputFile, Int_t Nevent = 5)
{
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
	// Vertex fit pulls
	Int_t Nbin = 100;
	TH1D* hXpull = new TH1D("hXpull", "Pull X vertex component", Nbin, -10., 10.);
	TH1D* hYpull = new TH1D("hYpull", "Pull Y vertex component", Nbin, -10., 10.);
	TH1D* hZpull = new TH1D("hZpull", "Pull Z vertex component", Nbin, -10., 10.);
	TH1D* hChi2 = new TH1D("hChi2", "Vertex #chi^{2}/N_{dof}", Nbin, 0., 10.);
	// Track vertex associations
	TH1D* hTrPrim = new TH1D("hTrPrim", "Available primary tracks", 41, -0.5, 40.5);
	TH1D* hTrFound = new TH1D("hTrFound", "Found primary tracks", 41, -0.5, 40.5);
	TH1D* hTrFoundG = new TH1D("hTrFoundG", "Found GOOD primary tracks", 41, -0.5, 40.5);
	TH1D* hTrFoundB = new TH1D("hTrFoundB", "Found BAD  primary tracks", 41, -0.5, 40.5);
	//
	// Loop over all events
	Int_t Nev = TMath::Min(Nevent, (Int_t)numberOfEntries);
	for (Int_t entry = 0; entry < Nev; ++entry)
	{
		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);
		Int_t Ntr = 0;	// # of starting tracks used in primary vertex
		Int_t NtrG = branchTrack->GetEntries();
		//std::cout << "Event opened containing " << NtrG << " tracks" << std::endl;
		TVectorD** pr = new TVectorD * [NtrG];
		TMatrixDSym** cv = new TMatrixDSym * [NtrG];
		//
		// test Particle branch
		Int_t Ngen = branchGenPart->GetEntries();
		//std::cout << "Nr. of generated particles: " << Ngen << std::endl;
		// If event contains at least 1 track
		//
		Double_t Nprim = 0.0;
		if (branchTrack->GetEntries() > 0)
		{
			// Loop on tracks
			for (Int_t it = 0; it < branchTrack->GetEntries(); it++)
			{
				Track* trk = (Track*)branchTrack->At(it);
				//
				// Start fitting all available tracks
				//
				// Reconstructed track parameters
				Double_t obsD0 = trk->D0;
				Double_t obsPhi = trk->Phi;
				Double_t obsC = trk->C;
				Double_t obsZ0 = trk->DZ;
				Double_t obsCtg = trk->CtgTheta;
				//std::cout << "Got track parameters for track " << it << std::endl;
				//
				// Load tracks for vertex fit if impact parameters is not ridiculous
				Double_t Dmax = 1.0;	// max is 1 mm
				if (TMath::Abs(obsD0) < Dmax) {
					Double_t oPar[5] = { obsD0, obsPhi, obsC, obsZ0, obsCtg };
					TVectorD obsPar(5, oPar);	// Fill observed parameters
					pr[Ntr] = new TVectorD(obsPar);
					cv[Ntr] = new TMatrixDSym(trk->CovarianceMatrix());
					Ntr++;
					//std::cout << "Track loaded Ntr= " << Ntr << std::endl;
				}
				//
				// Get associated generated particle
				GenParticle* gp = (GenParticle*)trk->Particle.GetObject();
				//std::cout << "GenParticle pointer "<<gp << std::endl;
				//
				// Position of origin in mm
				Double_t x = gp->X;
				Double_t y = gp->Y;
				Double_t z = gp->Z;
				//std::cout << "Got position of origin " << std::endl;
				//
				// Count tracks originating from the primary vertex
				if (x == 0.0 && y == 0.0) Nprim++;

			}		// End loop on tracks
		}
		//
		// Fit primary vertex
		//
		Int_t Nfound = Ntr;
		//std::cout << "Found tracks "<<Nfound << std::endl;
		Int_t MinTrk = 2;	// Minumum # tracks for vertex fit
		Double_t MaxChi2 = 6.;
		if (Ntr >= MinTrk) {
			VertexFit* Vtx = new VertexFit(Ntr, pr, cv);
			//std::cout << "Vertex fit created " << std::endl;
			//
			// Remove tracks with large chi2
			Bool_t Done = kFALSE;
			while (!Done){
				//std::cout << "After while " << std::endl;
				// Find largest Chi2 contribution
				TVectorD Chi2List = Vtx->GetVtxChi2List();	// Get contributions to Chi2
				//std::cout << "After Chi2List.  " << std::endl; Chi2List.Print();
				Double_t* Chi2L = new Double_t[Nfound];
				Chi2L = Chi2List.GetMatrixArray();
				Int_t iMax = TMath::LocMax(Nfound, Chi2L);
				//std::cout << "iMax =  "<<iMax << std::endl;
				Double_t Chi2Mx = Chi2L[iMax];
				//std::cout << "Chi2Mx "<<Chi2Mx << std::endl;
				if (Chi2Mx > MaxChi2 && Nfound > 2) {
					Vtx->RemoveTrk(iMax); 
					Nfound--;
				}
				else {
					Done = kTRUE;
				}
				//std::cout << "Done =  " << Done << std::endl;
				//delete[] Chi2L;
				//std::cout << "Array Chi2L removed " << std::endl;
			}
			//
			//std::cout << "Before getting vertex " << std::endl;
			//
			// Require minimum number of tracks in vertex
			Int_t Nmin = 4;
			if (Nfound >= Nmin) {
				TVectorD xvtx = Vtx->GetVtx();
				//std::cout << "Found vertex " << xvtx(0)<<", "<<xvtx(1)<<", "<<xvtx(2) << std::endl;
				TMatrixDSym covX = Vtx->GetVtxCov();
				Double_t Chi2 = Vtx->GetVtxChi2();
				Double_t Ndof = 2 * (Double_t)Nfound - 3;
				Double_t PullX = xvtx(0) / TMath::Sqrt(covX(0, 0));
				Double_t PullY = xvtx(1) / TMath::Sqrt(covX(1, 1));
				Double_t PullZ = xvtx(2) / TMath::Sqrt(covX(2, 2));
				//
				// Fill histograms
				hXpull->Fill(PullX);
				hYpull->Fill(PullY);
				hZpull->Fill(PullZ);
				hChi2->Fill(Chi2 / Ndof);
				//
				hTrPrim->Fill(Nprim);
				hTrFound->Fill((Double_t)Nfound);
				//std::cout << "Histograms filled " << std::endl;
			}
			//
			delete Vtx;
		}

		//std::cout << "Vertex chi2/Ndof = " << Chi2 / Ndof << std::endl;
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
	TCanvas* CnvN = new TCanvas("CnvN", "Primary tracks found", 100, 100, 900, 500);
	CnvN->Divide(2, 1);
	CnvN->cd(1);
	hTrPrim->Draw();
	CnvN-> cd(2);
	hTrFound->SetLineColor(kRed);
	hTrFound->Draw();
}
