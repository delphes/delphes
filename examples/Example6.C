/*
Simple macro showing how to access branches from the delphes output root file,
loop over events, and plot track pulls and geometrical acceptance.

root -l examples/Example6.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "modules/TrackCovariance.h"
#include "external/TrackCovariance/TrkUtil.h"
#endif


//------------------------------------------------------------------------------

void Example6(const char* inputFile)
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

	// Book histograms
	//
	// Generated track parameters
	TH1* histDgen = new TH1F("h_Dgen", "Generated impact parameter", 100, -2., 2.);			// mm
	TH1* histP0gen = new TH1F("h_P0gen", "Generated #phi_{0}", 100, -TMath::Pi(), TMath::Pi());	// rad
	TH1* histCgen = new TH1F("h_Cgen", "Generated half curvature", 200, -.01, .01);				// mm
	TH1* histZ0gen = new TH1F("h_Z0gen", "Generated Z_{0}", 100, -2., 2.);							// mm
	TH1* histCtgen = new TH1F("h_Ctgen", "Generated cotg(#theta)", 100, -10., 10.);
	//
	// Reconstructed impact parameters
	TH1* histDobs = new TH1F("h_Dobs", "Reconstructed impact parameter", 100, -2., 2.);			// mm
	TH1* histP0obs = new TH1F("h_P0obs", "Reconstructed #phi_{0}", 100, -TMath::Pi(), TMath::Pi());	// rad
	TH1* histCobs = new TH1F("h_Cobs", "Reconstructed half curvature", 200, -.01, .01);				// mm
	TH1* histZ0obs = new TH1F("h_Z0obs", "Reconstructed Z_{0}", 100, -2., 2.);							// mm
	TH1* histCtobs = new TH1F("h_Ctobs", "Reconstructed cotg(#theta)", 100, -10., 10.);
	//
	// Track parameter pulls
	TH1* histDpull = new TH1F("h_Dpull", "Pull impact parameter", 100, -10., 10.);
	TH1* histP0pull = new TH1F("h_P0pull", "Pull #phi_{0}", 100, -10., 10.);
	TH1* histCpull = new TH1F("h_Cpull", "Pull half curvature", 100, -10., 10.);
	TH1* histZ0pull = new TH1F("h_Z0pull", "Pull Z_{0}", 100, -10., 10.);
	TH1* histCtpull = new TH1F("h_Ctpull", "Pull cotg(#theta)", 100, -10., 10.);
	//
	// Accptance plots
	Double_t minPtAcc = 1.0;	// 1 GeV
	Double_t maxCtAcc = 1.0;	// 45 degrees
	// All tracks
	TH1* histPtgen = new TH1F("h_Ptgen", "Generated Pt", 500, 0., 50.);
	TH1* histPtobs = new TH1F("h_Ptobs", "Reconstructed Pt", 500, 0., 50.);
	// CEntral track over min Pt
	TH1* histPtgenC = new TH1F("h_PtgenC", "Generated Pt - Central", 500, 0., 50.);		// pt for central tracks;
	TH1* histPtobsC = new TH1F("h_PtobsC", "Reconstructed Pt - Central", 500, 0., 50.);	// pt for central
	TH1* histCtgenH = new TH1F("h_CtgenH", "Generateded Cotangent", 100, -10., 10.);		// cot(theta) for high pt tracks;
	TH1* histCtobsH = new TH1F("h_CtobsH", "Reconstructed Cotangent", 100, -10., 10.);	// cot(theta) for high pt tracks
	// All tracks
	TH1* histAccCtg = new TH1F("h_AccCtg", "Cotangent acceptance", 100, -10., 10.);
	TH1* histAccPt = new TH1F("h_AccPt", "Pt acceptance", 500, 0., 50.);
	// Tracks over min pt
	TH1* histAccCtgH = new TH1F("h_AccCtgH", "Cotangent acceptance (high pt)", 100, -10., 10.);
	// tracks above min angle
	TH1* histAccPtC = new TH1F("h_AccPtC", "Pt acceptance (central tracks)", 500, 0., 50.);

	//
	// Get magnetifc field
	Double_t Bz = treeReader->GetInfo("Bz");

	// Loop over all events
	for (Int_t entry = 0; entry < numberOfEntries; ++entry)
	{
		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);

		// If event contains at least 1 track
		//
		if (branchTrack->GetEntries() > 0)
		{
			// Loop on tracks
			for (Int_t it = 0; it < branchTrack->GetEntries(); it++)
			{
				Track* trk = (Track*)branchTrack->At(it);
				//
				// Reconstructed track parameters
				Double_t obsD0  = trk->D0;
				Double_t obsPhi = trk->Phi;
				Double_t obsC   = trk->C;
				Double_t obsZ0  = trk->DZ;
				Double_t obsCtg = trk->CtgTheta;
				// Fill histograms
				histDobs ->Fill(obsD0);
				histP0obs->Fill(obsPhi);
				histCobs ->Fill(obsC);
				histZ0obs->Fill(obsZ0);
				histCtobs->Fill(obsCtg);
				//
				// Get associated generated particle
				GenParticle* gp = (GenParticle*)trk->Particle.GetObject();
				//
				// Position of origin in meters
				Double_t x = 1.0e-3 * gp->X;
				Double_t y = 1.0e-3 * gp->Y;
				Double_t z = 1.0e-3 * gp->Z;
				TVector3 xv(x, y, z);
				//
				// Momentum at origin
				Double_t px = gp->Px;
				Double_t py = gp->Py;
				Double_t pz = gp->Pz;
				TVector3 pv(px, py, pz);
				//
				// Unsmeared original parameters
				Double_t Q = gp->Charge;
				TVectorD gPar   = TrkUtil::XPtoPar(xv, pv, Q, Bz);	// Get parameters from position and momentum
				TVectorD genPar = TrkUtil::ParToMm(gPar);			// Convert to mm
				// Unpack parameters
				Double_t genD0  = genPar(0);
				Double_t genPhi = genPar(1);
				Double_t genC   = genPar(2);
				Double_t genZ0  = genPar(3);
				Double_t genCtg = genPar(4);
				Double_t genPt  = pv.Pt();
				//
				// Calculate and fill pulls
				Double_t pullD0  = (obsD0  - genD0)  / trk->ErrorD0;
				Double_t pullPhi = (obsPhi - genPhi) / trk->ErrorPhi;
				Double_t pullC   = (obsC   - genC)   / trk->ErrorC;
				Double_t pullZ0  = (obsZ0  - genZ0)  / trk->ErrorDZ;
				Double_t pullCtg = (obsCtg - genCtg) / trk->ErrorCtgTheta;
				//
				histDpull ->Fill(pullD0);
				histP0pull->Fill(pullPhi);
				histCpull ->Fill(pullC);
				histZ0pull->Fill(pullZ0);
				histCtpull->Fill(pullCtg);
				//
				// Acceptance plots
				Double_t obsPt = trk->PT;
				histPtobs->Fill(obsPt);
				if (TMath::Abs(genCtg) < maxCtAcc)histPtobsC->Fill(obsPt); // pt for central tracks
				if (obsPt > minPtAcc)histCtobsH->Fill(obsCtg);				// cot(theta) for high pt tracks
			}		// End loop on tracks

			// If event contains at least 1 generated charged
			//
			if (branchGenPart->GetEntries() > 0)
			{
				// Loop on generated particles
				for (Int_t it = 0; it < branchGenPart->GetEntries(); it++) {
					GenParticle* gpart = (GenParticle*)branchGenPart->At(it);
					//
					// Plot charged particle parameters
					// Only final state particles (Status = 1)
					if (gpart->Status == 1 && TMath::Abs(gpart->Charge) > 0) {
						//
						// Position of origin in meters
						Double_t x = 1.0e-3 * gpart->X;
						Double_t y = 1.0e-3 * gpart->Y;
						Double_t z = 1.0e-3 * gpart->Z;
						TVector3 xv(x, y, z);
						//
						// Momentum at origin
						Double_t px = gpart->Px;
						Double_t py = gpart->Py;
						Double_t pz = gpart->Pz;
						TVector3 pv(px, py, pz);
						//
						// Original parameters
						Double_t Q = gpart->Charge;
						TVectorD gPar = TrkUtil::XPtoPar(xv, pv, Q, Bz);	// Get parameters from position and momentum
						TVectorD genPar = TrkUtil::ParToMm(gPar);			// Convert to mm
						// Unpack parameters
						Double_t genD0  = genPar(0);
						Double_t genPhi = genPar(1);
						Double_t genC   = genPar(2);
						Double_t genZ0  = genPar(3);
						Double_t genCtg = genPar(4);
						Double_t genPt  = pv.Pt();
						// Fill histograms
						histDgen ->Fill(genD0);
						histP0gen->Fill(genPhi);
						histCgen ->Fill(genC);
						histZ0gen->Fill(genZ0);
						histCtgen->Fill(genCtg);
						//
						// Acceptance plots
						histPtgen->Fill(genPt);
						if (TMath::Abs(genCtg) < maxCtAcc)histPtgenC->Fill(genPt);
						if (genPt > minPtAcc)histCtgenH->Fill(genCtg);
					}
				}
			}
		}
	}
	//
	// Show resulting histograms
	//
	TCanvas* Cnv = new TCanvas("Cnv", "Delphes generated track plots", 50, 50, 900, 500);
	Cnv->Divide(3, 2);
	Cnv->cd(1); histDgen->Draw();
	Cnv->cd(2); histP0gen->Draw();
	Cnv->cd(3); histCgen->Draw();
	Cnv->cd(4); histZ0gen->Draw();
	Cnv->cd(5); histCtgen->Draw();
	//
	TCanvas* Cnv1 = new TCanvas("Cnv1", "Delphes observed track plots", 100, 100, 900, 500);
	Cnv1->Divide(3, 2);
	Cnv1->cd(1); histDobs->Draw();
	Cnv1->cd(2); histP0obs->Draw();
	Cnv1->cd(3); histCobs->Draw();
	Cnv1->cd(4); histZ0obs->Draw();
	Cnv1->cd(5); histCtobs->Draw();
	//
	TCanvas* Cnv2 = new TCanvas("Cnv2", "Delphes observed track pulls", 150, 150, 900, 500);
	Cnv2->Divide(3, 2);
	Cnv2->cd(1); gPad->SetLogy(1);
	histDpull->Draw();
	Cnv2->cd(2); gPad->SetLogy(1);
	histP0pull->Draw();
	Cnv2->cd(3); gPad->SetLogy(1);
	histCpull->Draw();
	Cnv2->cd(4); gPad->SetLogy(1);
	histZ0pull->Draw();
	Cnv2->cd(5); gPad->SetLogy(1);
	histCtpull->Draw();
	//
	// Acceptance plots
	//
	TCanvas* Cnv4 = new TCanvas("Cnv4", "Delphes acceptance", 200, 200, 500, 500);
	Cnv4->Divide(2, 2);
	Cnv4->cd(1);
	histCtgen->Draw();
	histCtobs->SetLineColor(kRed);
	histCtobs->Draw("SAME");
	Cnv4->cd(2);
	Int_t NbCt = histCtgen->GetNbinsX();
	for (Int_t i = 1; i < NbCt + 1; i++)
	{
		Float_t cgen = histCtgen->GetBinContent(i);
		Float_t cobs = histCtobs->GetBinContent(i);
		Float_t AccCtg = 0.0;
		if (cgen > 0.0) AccCtg = cobs / cgen;
		histAccCtg->SetBinContent(i, AccCtg);
	}
	histAccCtg->Draw();
	Cnv4->cd(3); gPad->SetLogy(1);
	histPtgen->Draw();
	histPtobs->SetLineColor(kRed);
	histPtobs->Draw("SAME");
	Cnv4->cd(4);
	Int_t NbPt = histPtgen->GetNbinsX();
	for (Int_t i = 1; i < NbPt + 1; i++)
	{
		Float_t pgen = histPtgen->GetBinContent(i);
		Float_t pobs = histPtobs->GetBinContent(i);
		Float_t AccPt = 0.0;
		if (pgen > 0.0) AccPt = pobs / pgen;
		histAccPt->SetBinContent(i, AccPt);
	}
	histAccPt->Draw();
	//
	TCanvas* Cnv5 = new TCanvas("Cnv5", "Delphes acceptance (constrained)", 250, 250, 500, 500);
	Cnv5->Divide(2, 2);
	Cnv5->cd(1);
	histCtgenH->Draw();
	histCtobsH->SetLineColor(kRed);
	histCtobsH->Draw("SAME");
	Cnv5->cd(2);
	NbCt = histCtgenH->GetNbinsX();
	for (Int_t i = 1; i < NbCt + 1; i++)
	{
		Float_t cgen = histCtgenH->GetBinContent(i);
		Float_t cobs = histCtobsH->GetBinContent(i);
		Float_t AccCtg = 0.0;
		if (cgen > 0.0) AccCtg = cobs / cgen;
		histAccCtgH->SetBinContent(i, AccCtg);
	}
	histAccCtgH->Draw();
	Cnv5->cd(3); gPad->SetLogy(1);
	histPtgenC->Draw();
	histPtobsC->SetLineColor(kRed);
	histPtobsC->Draw("SAME");
	Cnv5->cd(4);
	NbPt = histPtgenC->GetNbinsX();
	for (Int_t i = 1; i < NbPt + 1; i++)
	{
		Float_t pgen = histPtgenC->GetBinContent(i);
		Float_t pobs = histPtobsC->GetBinContent(i);
		Float_t AccPt = 0.0;
		if (pgen > 0.0) AccPt = pobs / pgen;
		histAccPtC->SetBinContent(i, AccPt);
	}
	histAccPtC->Draw();
}
