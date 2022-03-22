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
#include <TVectorD.h>
#include <TMath.h>

#endif



//------------------------------------------------------------------------------

void ExampleKsVtxFit(const char* inputFile, Int_t Nevent = 5)
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
	Int_t Nbin = 100;
	TH1D* hXpull = new TH1D("hXpull", "Pull X vertex component", Nbin, -10., 10.);
	TH1D* hYpull = new TH1D("hYpull", "Pull Y vertex component", Nbin, -10., 10.);
	TH1D* hZpull = new TH1D("hZpull", "Pull Z vertex component", Nbin, -10., 10.);
	TH1D* hChi2 = new TH1D("hChi2", "Vertex #chi^{2}/N_{dof}", Nbin, 0., 10.);
	//
	// Loop over all events
	Int_t Nev = TMath::Min(Nevent, (Int_t)numberOfEntries);
	for (Int_t entry = 0; entry < Nev; ++entry)
	{
		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);
		// 
		Int_t NtrG = branchTrack->GetEntries();
		const Int_t NtFit = 2;
		TVectorD** pr = new TVectorD * [NtFit];
		TMatrixDSym** cv = new TMatrixDSym * [NtFit];
		std::cout<<"Start of event "<<entry<<std::endl;
		//
		// Search for generated Ks--> pi+ pi-
		Int_t KsID = 310;
		Int_t NKs = 0;
		std::vector<GenParticle*> gpi1;
		std::vector<GenParticle*> gpi2;
		std::vector<Track*> gti1;
		std::vector<Track*> gti2;
		std::vector<TVectorD*> gKsVtx;
		std::vector<TVectorD*> rKsVtx;
		//
		if (branchGenPart->GetEntries() > 1)
		{
		    for(Int_t ig=0; ig<branchGenPart->GetEntries(); ig++)
		    {
			GenParticle* gp = (GenParticle*)branchGenPart->At(ig);
			if(gp->PID == KsID)	// Found Ks
			{
				// Store daughters
				gpi1.push_back( (GenParticle*)branchGenPart->At(gp->D1));
				gpi2.push_back( (GenParticle*)branchGenPart->At(gp->D2));
				Track* Tnull = 0;
				gti1.push_back(Tnull);
				gti2.push_back(Tnull);
				// Store vertex
				Double_t Vnull[3] = {0.,0.,0.};
				TVectorD rVnull(3,Vnull);
				Double_t Vgen[3] = {gpi1[NKs]->X,gpi1[NKs]->Y,gpi1[NKs]->Z};
				TVectorD gVert(3,Vgen);
				gKsVtx.push_back(new TVectorD(gVert));
				rKsVtx.push_back(new TVectorD(rVnull));
				//
				NKs++;
			}
		    }
		}	// Done searching
		std::cout<<"Done searching. Found "<<NKs<<" K0s"<<std::endl;
		//
		// If event contains at least 1 Ks and 1 track
		//
		if (branchTrack->GetEntries() > 1 && NKs > 0)
		{
			// Loop on tracks
			for (Int_t it = 0; it < branchTrack->GetEntries(); it++)
			{
				//std::cout<<"Start track loop. it ="<<it<<std::endl;
				Track* trk = (Track*)branchTrack->At(it);
				//
				// Get associated generated particle
				GenParticle* gt = (GenParticle*)trk->Particle.GetObject();
				Int_t mp = gt->M1;	// Mother
				GenParticle* gm = (GenParticle*)branchGenPart->At(mp);
				if(gm->PID == KsID)	// Mother is Ks
				{
				   for(Int_t k=0; k<NKs; k++)
				   {
					if(gt == gpi1[k]) gti1[k] = trk;
					if(gt == gpi2[k]) gti2[k] = trk;
				   }
				}
				//std::cout<<"End track loop"<<std::endl;
			}		// End loop on tracks
			//
			//	Look for K0s to fit
			//
			for(Int_t k=0; k<NKs; k++)	//Scan Ks
			{
				Bool_t found = (gti1[k] != 0) && (gti2[k] != 0);
				//std::cout<<"Prepare Ks # "<<k<<" for fitting. Found = "<<found<<std::endl;
				TVectorD gVtx = *gKsVtx[k];
				if(gti1[k] != 0 && gti2[k] != 0)	// if both tracks found
				{	
					//std::cout<<"Both tracks found"<<std::endl;
					// First leg
					Double_t obsD0 = gti1[k]->D0;
					Double_t obsPhi = gti1[k]->Phi;
					Double_t obsC = gti1[k]->C;
					Double_t obsZ0 = gti1[k]->DZ;
					Double_t obsCtg = gti1[k]->CtgTheta;
					// 
					Double_t oPar1[5] = { obsD0, obsPhi, obsC, obsZ0, obsCtg };
					TVectorD obsPar1(5, oPar1);	// Fill observed parameters
					pr[0] = new TVectorD(obsPar1);
					cv[0] = new TMatrixDSym(gti1[k]->CovarianceMatrix());
					// Second leg
					obsD0 = gti2[k]->D0;
					obsPhi = gti2[k]->Phi;
					obsC = gti2[k]->C;
					obsZ0 = gti2[k]->DZ;
					obsCtg = gti2[k]->CtgTheta;
					// 
					Double_t oPar2[5] = { obsD0, obsPhi, obsC, obsZ0, obsCtg };
					TVectorD obsPar2(5, oPar2);	// Fill observed parameters
					pr[1] = new TVectorD(obsPar2);
					cv[1] = new TMatrixDSym(gti2[k]->CovarianceMatrix());
					//
					// Fit vertex
					//
					Double_t Rmin = 17.;	// Lowest layer
					Double_t Rvtx = TMath::Sqrt(gVtx(0)*gVtx(0)+gVtx(0)*gVtx(0));
					VertexFit* Vtx = new VertexFit(NtFit, pr, cv);
					if(Rvtx > Rmin)Vtx->SetStartR(Rvtx+10.);
					//std::cout<<"BEFORE VERTEX FIT. Rvtx = "<<Rvtx<<std::endl;
					TVectorD xvtx = Vtx->GetVtx();
					TMatrixDSym covX = Vtx->GetVtxCov();
					Double_t Chi2 = Vtx->GetVtxChi2();
					//
					// Fill plots
					//
					Double_t PullX = (xvtx(0)-gVtx(0)) / TMath::Sqrt(covX(0, 0));
					Double_t PullY = (xvtx(1)-gVtx(1)) / TMath::Sqrt(covX(1, 1));
					Double_t PullZ = (xvtx(2)-gVtx(2)) / TMath::Sqrt(covX(2, 2));
					hXpull->Fill(PullX);
					hYpull->Fill(PullY);
					hZpull->Fill(PullZ);
					hChi2->Fill(Chi2);
					std::cout<<"xg: "<<gVtx(0)<<", yg: "<<gVtx(1)<<", zg: "<<gVtx(2)<<std::endl;
					std::cout<<"xr: "<<xvtx(0)<<", yr: "<<xvtx(1)<<", zr: "<<xvtx(2)<<std::endl;
					// Cleanup
					for (Int_t i = 0; i < NtFit; i++) delete pr[i];
					for (Int_t i = 0; i < NtFit; i++) delete cv[i];
				}
			} // End loop on Ks
			//std::cout<<"Joust out of Ks loop"<<std::endl;
			//
			// Clean
			delete[] pr;
			delete[] cv;
			for(Int_t k=0; k<NKs; k++)
			{
				delete gKsVtx[k];
				delete rKsVtx[k];
			}
			gpi1.clear();;
			gpi2.clear();
			gti1.clear();
			gti2.clear();
			gKsVtx.clear();
			rKsVtx.clear();
			//std::cout<<"Finished cleanup"<<std::endl;
		} // End if on Ks present
	}
	//
	// Show resulting histograms
	//
	TCanvas* Cnv = new TCanvas("Cnv", "Delphes Ks vertex pulls", 50, 50, 900, 500);
	Cnv->Divide(2, 2);
	Cnv->cd(1); gPad->SetLogy(1); gStyle->SetOptFit(1111);
	hXpull->Fit("gaus"); hXpull->Draw();
	Cnv->cd(2); gPad->SetLogy(1); gStyle->SetOptFit(1111);
	hYpull->Fit("gaus"); hYpull->Draw();
	Cnv->cd(3); gPad->SetLogy(1); gStyle->SetOptFit(1111);
	hZpull->Fit("gaus"); hZpull->Draw();
	Cnv->cd(4); hChi2->Draw();
}
