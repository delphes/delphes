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
	TH1D* hTrDiff = new TH1D("hTrDiff", "Found - available tracks", 21, -10.5, 10.5);
	//
	// Loop over all events
	Int_t Nev = TMath::Min(Nevent, (Int_t)numberOfEntries);
	for (Int_t entry = 0; entry < Nev; ++entry)
	{
		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);
		Int_t NtrG = branchTrack->GetEntries();
		if(entry%500 ==0)std::cout << "Event "<<entry<<" opened containing " << NtrG << " tracks" << std::endl;
		TVectorD** pr = new TVectorD * [NtrG];
		TMatrixDSym** cv = new TMatrixDSym * [NtrG];
		//
		// test Particle branch
		Int_t Ngen = branchGenPart->GetEntries();
		//std::cout << "Nr. of generated particles: " << Ngen << std::endl;
		// If event contains at least 1 track
		//
		Double_t Nprim = 0.0;
		Double_t xpv = 0.0;		// Init true primary vertex position
		Double_t ypv = 0.0;
		Double_t zpv = 0.0;
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
				// Load all tracks for vertex fit 
				Double_t oPar[5] = { obsD0, obsPhi, obsC, obsZ0, obsCtg };
				TVectorD obsPar(5, oPar);	// Fill observed parameters
				pr[it] = new TVectorD(obsPar);
				cv[it] = new TMatrixDSym(trk->CovarianceMatrix());
				// Check covariance
				/*
				TMatrixDSymEigen Eign(*cv[it]);
				TVectorD lambda = Eign.GetEigenValues();
				Bool_t good = kTRUE;
				for(Int_t i=0; i<5; i++) if(lambda(i)<0.0) good = kFALSE;
				if(!good){
					std::cout<<"Cov["<<it<<"]"<<" eigenvalues: "
					<<lambda(0)<<", "<<lambda(1)<<", "<<lambda(2)<<", "
					<<lambda(3)<<", "<<lambda(4)<<std::endl;
					TMatrixDSym Ctmp = *cv[it];
					std::cout<<"Error D (dir - cov) "<<trk->ErrorD0
					<<" - "<<TMath::Sqrt(Ctmp(0,0))<<std::endl;
				}
				*/
				//std::cout << "Track loaded it= " << it << std::endl;
			//
			// Find true primary vertex
				GenParticle* gp = (GenParticle*)trk->Particle.GetObject();
				//std::cout << "GenParticle pointer "<<gp << std::endl;
				//
				// Position of origin in mm
				Double_t x = gp->X;
				Double_t y = gp->Y;
				Double_t z = gp->Z;
				Bool_t prim = kTRUE;	// Is primary?
				Int_t mp = gp->M1;	// Mother
				while (mp > 0) {
					GenParticle* gm =
						(GenParticle*)branchGenPart->At(mp);
					Double_t xm = gm->X;
					Double_t ym = gm->Y;
					Double_t zm = gm->Z;
					if (x != xm || y != ym || z != zm) {
						prim = kFALSE;
						break;
					}
					else mp = gm->M1;
				}
				if (prim) {		// It's a primary track
					Nprim++;
					xpv = x;	// Store true primary
					ypv = y;
					zpv = z;
				}

			}		// End loop on tracks
		}
		std::cout<<"PVtxFind true vertex: Nprim= "<<Nprim<<", x,y,z= "<<xpv<<", "<<ypv<<", "<<zpv<<std::endl;
		//
		// Find primary vertex
		//
		//Beam constraint
		TVectorD xpvc(3);
		xpvc(0) = 1.0;
		xpvc(1) = -2.0;
		xpvc(2) = 10.0;
		TMatrixDSym covpvc(3); covpvc.Zero();
		covpvc(0, 0) = 0.0097 * 0.0097;
		covpvc(1, 1) = 2.55e-05 * 2.55e-05;
		covpvc(2, 2) = 0.64 * 0.64;
		if(Nprim == 0){
			xpv = xpvc(0);
			ypv = xpvc(1);
			zpv = xpvc(2);
		}
		//
		//
		// Skim tracks
		Int_t nSkim = 0;
		Int_t* nSkimmed = new Int_t[NtrG];
		TVectorD** PrSk = new TVectorD * [1];
		TMatrixDSym** CvSk = new TMatrixDSym * [1];
		Double_t MaxChi2 = 2.0+3*2.;
		for (Int_t n = 0; n < NtrG; n++) {
			PrSk[0] = new TVectorD(*pr[n]);
			CvSk[0] = new TMatrixDSym(*cv[n]);
			VertexFit* Vskim = new VertexFit(1,PrSk, CvSk);
			Vskim->AddVtxConstraint(xpvc, covpvc);
			Double_t Chi2One = Vskim->GetVtxChi2();
			//std::cout<<"Track "<<n<<", Chi2 = "<<Chi2One<<std::endl;
			if (Chi2One < MaxChi2) {
				nSkimmed[nSkim] = n;
				//std::cout << "nSkimmed[" << nSkim << "] = " << n << std::endl;
				nSkim++;
			}
			// Cleanup
			delete Vskim;
		}
		delete PrSk[0];
		delete CvSk[0];
		delete[] PrSk;
		delete[] CvSk;
		//
		// Load tracks for primary fit
		Int_t MinTrk = 1;	// Minumum # tracks for vertex fit
		if (nSkim >= MinTrk) {
			TVectorD** PrFit = new TVectorD * [nSkim];
			TMatrixDSym** CvFit = new TMatrixDSym * [nSkim];
			for (Int_t n = 0; n < nSkim; n++) {
				PrFit[n] = new TVectorD(*pr[nSkimmed[n]]);
				CvFit[n] = new TMatrixDSym(*cv[nSkimmed[n]]);
			}
			delete[] nSkimmed;
			Int_t Nfound = nSkim;
			if(entry%500 ==0)std::cout << "Found tracks "<<Nfound << std::endl;
			const Int_t MaxFound = 100; Double_t Chi2LL[MaxFound]; Double_t *Chi2L = Chi2LL;
			VertexFit* Vtx = new VertexFit(nSkim, PrFit, CvFit);
			//std::cout << "Vertex fit created " << std::endl;
			Vtx->AddVtxConstraint(xpvc, covpvc);
			//
			// Remove tracks with large chi2
			Double_t MaxChi2Fit = 4.5;
			Bool_t Done = kFALSE;
			while (!Done) {
				//std::cout << "After while " << std::endl;
				// Find largest Chi2 contribution
				TVectorD Chi2List = Vtx->GetVtxChi2List();	// Get contributions to Chi2
				//std::cout << "After Chi2List.  " << std::endl; Chi2List.Print();
				//Double_t* Chi2L = new Double_t[Nfound];
				Chi2L = Chi2List.GetMatrixArray();
				Int_t iMax = TMath::LocMax(Nfound, Chi2L);
				//std::cout << "iMax =  "<<iMax << std::endl;
				Double_t Chi2Mx = Chi2L[iMax];
				//std::cout << "Chi2Mx "<<Chi2Mx << std::endl;
				if (Chi2Mx > MaxChi2Fit && Nfound > 1) {
					//std::cout << "Before remove.  Nfound = "<<Nfound << std::endl;
					Vtx->RemoveTrk(iMax);
					//std::cout << "After remove." << std::endl;
					Nfound--;
				}
				else {
					Done = kTRUE;
				}
			}
			//
			//std::cout << "Before getting vertex " << std::endl;
			//
			// Require minimum number of tracks in vertex
			Int_t Nmin = 1;
			if (Nfound >= Nmin) {
				TVectorD xvtx = Vtx->GetVtx();
				std::cout << "Found vertex " << xvtx(0)<<", "<<xvtx(1)<<", "<<xvtx(2) << std::endl;
				TMatrixDSym covX = Vtx->GetVtxCov();
				Double_t Chi2 = Vtx->GetVtxChi2();
				Double_t Ndof = 2 * (Double_t)Nfound;
				Double_t PullX = (xvtx(0) - xpv) / TMath::Sqrt(covX(0, 0));
				Double_t PullY = (xvtx(1) - ypv) / TMath::Sqrt(covX(1, 1));
				Double_t PullZ = (xvtx(2) - zpv) / TMath::Sqrt(covX(2, 2));
				//
				// Fill histograms
				hXpull->Fill(PullX);
				hYpull->Fill(PullY);
				hZpull->Fill(PullZ);
				hChi2->Fill(Chi2 / Ndof);
				//
				hTrPrim->Fill(Nprim);
				hTrFound->Fill((Double_t)Nfound);
				hTrDiff->Fill((Double_t)Nfound-Nprim);
				//std::cout << "Histograms filled " << std::endl;
			}
			//
			// Clean
			delete Vtx;
			for (Int_t i = 0; i < nSkim; i++) delete PrFit[i];
			for (Int_t i = 0; i < nSkim; i++) delete CvFit[i];
			delete[] PrFit;
			delete[] CvFit;
		}

		//std::cout << "Vertex chi2/Ndof = " << Chi2 / Ndof << std::endl;
		//
		// Cleanup
		for (Int_t i = 0; i < NtrG; i++) delete pr[i];
		for (Int_t i = 0; i < NtrG; i++) delete cv[i];
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
	CnvN->Divide(2, 2);
	CnvN->cd(1);
	hTrPrim->Draw();
	CnvN->cd(2);
	hTrFound->SetLineColor(kRed);
	hTrFound->Draw();
	CnvN->cd(3);
	hTrDiff->Draw();
}
