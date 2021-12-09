#include <iostream>
#include <TMath.h>
#include <TVectorD.h>
#include <TVector3.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraph.h>
#include <TLine.h>
#include <TRandom.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include "external/TrackCovariance/TrkUtil.h"


void ExampleClCount(Int_t Nev = 1000, Double_t ang = 90., Int_t Opt=0)
{
	//
	// Constants
	Double_t angRad = ang*TMath::Pi() / 180.;	// To radians
	Double_t KpMass = 0.493677;		// Charged Kaon mass
	Double_t PiMass = 0.13957018;	// Charged Pion mass
	Double_t Bz = 2.0;				// Field
	TVector3 x(0.0, 0.0, 0.0);		// Track origin
	//
	Double_t pmin = 0.1;				// Momentum range
	Double_t pmax = 50.;
	//
	// Setup chamber
	//
	TrkUtil *TU = new TrkUtil(Bz);
	//cout << "TU done" << endl;
	Double_t Rmin = 0.35;
	Double_t Rmax = 2.0;
	Double_t Zmin = -2.0;
	Double_t Zmax = 2.0;
	TU->SetDchBoundaries(Rmin, Rmax, Zmin, Zmax);
	TU->SetGasMix(Opt);	// gas selection
	char title[50];
	Double_t y_min = 0.0;
	Double_t y_max = 0.0;
	switch (Opt)
	{
	case 0:
	{
		sprintf(title, "Helium 90 - Isobutane 10");	// He-Isobutane
		y_min = 1500.; y_max = 3500;
		break;
	}
	case 1:
	{
		sprintf(title, "Helium 100");				// pure He
		y_min = 400.; y_max = 1000;
		break;
	}
	case 2:
	{
		sprintf(title, "Argon 50 - Ethane 50");		// Argon - Ethane
		y_min = 6000.; y_max = 8500.;
		break;
	}
	case 3:
	{
		sprintf(title, "Argon 100");				// pure Argon
		y_min = 4000.; y_max = 6000.;
		break;
	}
	}
	//cout << "DchBoundaries set" << endl;
	//
	// Histograms
	//
	TH2F *h_kaon = new TH2F("h_kaon", "Nclusters distribution", 200, 0., 50., 50, y_min, y_max);
	TH2F *h_pion = new TH2F("h_pion", "Nclusters distribution", 200, 0., 50., 50, y_min, y_max);
	//cout << "histograms booked" << endl;
	//
	// Fill plots
	//
	for (Int_t n = 0; n < Nev; n++)
	{
		Double_t R = gRandom->Rndm();
		Double_t pval = pmin + R*(pmax - pmin);
		Double_t px = pval*TMath::Sin(angRad);
		Double_t pz = pval*TMath::Cos(angRad);
		TVector3 p(px, 0.0, pz);
		Double_t Q = 1.0;
		TVectorD Par = TrkUtil::XPtoPar(x, p, Q, Bz);
		//cout << "Parameters done. pval = "<< pval << endl;
		//
		Double_t Ncl = 0.0;
		if (TU->IonClusters(Ncl, KpMass, Par))
		{
			//cout << "Kaon extracted. pval = "<< pval << endl;
			h_kaon->Fill(pval, Ncl);
			//cout << "Kaon filled. Ncl = "<<Ncl << endl;
		}
		if (TU->IonClusters(Ncl, PiMass, Par))
		{
			//cout << "Pion extracted. pval = "<<pval << endl;
			h_pion->Fill(pval, Ncl);
			//cout << "Pion filled. Ncl = "<<Ncl << endl;
		}
	}
	//
	// graph K/pi separation
	//
	const Int_t Nint = 501;
	Double_t pmom[Nint];
	Double_t pmn = 0.2; Double_t pmx = 100.;
	Double_t stp = (TMath::Log(pmx) - TMath::Log(pmn)) / (Double_t)(Nint - 1);
	for (Int_t i = 0; i < Nint; i++)pmom[i] = TMath::Exp(i * stp + TMath::Log(pmn));
	Double_t SigDiff[Nint];
	for (Int_t i = 0; i < Nint; i++)
	{
		Double_t bgK = pmom[i] / KpMass;
		Double_t bgP = pmom[i] / PiMass;
		Double_t px = pmom[i]*TMath::Sin(angRad);
		Double_t pz = pmom[i]*TMath::Cos(angRad);
		TVector3 p(px, 0.0, pz);
		Double_t Q = 1.0;
		TVectorD Par = TrkUtil::XPtoPar(x, p, Q, Bz);
		Double_t tLen = TU->TrkLen(Par);
		Double_t CluKp = TrkUtil::Nclusters(bgK,Opt)*tLen;
		Double_t CluKpErr = TMath::Sqrt(CluKp);
		Double_t CluPi = TrkUtil::Nclusters(bgP,Opt)*tLen;
		Double_t CluPiErr = TMath::Sqrt(CluPi);
		//cout << "Momentum = " << pmom[i] << ", length= " << tLen
		//	<< ", CluKp = " << CluKp << ", CluPi = " << CluPi << endl;
		//if (CluKp + CluPi > 0.0)SigDiff[i] = TMath::Abs(CluPi - CluKp) / TMath::Sqrt(CluKp + CluPi);
		if (CluKp + CluPi > 0.0)SigDiff[i] = TMath::Abs(CluPi - CluKp) / (0.5*(CluKpErr + CluPiErr));
	}
	//
	// Plots
	//
	TCanvas *cnv = new TCanvas("cnv", title, 50, 50, 800, 500);
	cnv->Divide(2, 1);
	cnv->cd(1);
	gStyle->SetOptStat(0);
	h_kaon->SetMarkerColor(kRed);
	h_kaon->SetLineColor(kRed);
	h_kaon->GetXaxis()->SetTitle("Momentum (GeV)");
	h_kaon->Draw();
	h_pion->SetMarkerColor(kBlack);
	h_pion->SetLineColor(kBlack);
	h_pion->Draw("SAME");
	//
	TLegend* lg = new TLegend(0.1, 0.9, 0.3, 0.70);
	TString LgTitle = "Particle type:";
	lg->SetHeader(LgTitle);
	lg->AddEntry(h_pion, "#pi", "L");
	lg->AddEntry(h_kaon, "#color[2]{K}", "L");
	lg->Draw();
	cnv->cd(2);
	gPad->SetLogx();
	TGraph *gr = new TGraph(Nint, pmom, SigDiff);
	gr->SetMinimum(0.0);
	gr->SetMaximum(8.0);
	gr->GetXaxis()->SetLimits(0.,100.);
	gr->SetTitle("K/#pi separation in nr. of #sigma");
	gr->GetXaxis()->SetTitle("Momentum (GeV)");
	gr->SetLineColor(kBlue);
	gr->Draw("APC");
	TLine *line = new TLine(0.0, 3.0, 100., 3.0);
	line->Draw("SAME");
	//
	TCanvas* cnvf = new TCanvas("cnvf", "Interpolating function", 100, 100, 800, 500);
	cnvf->Divide(1, 1);
	cnvf->cd(1);
	gPad->SetLogx();
	TF1* f_ncl = new TF1("f_ncl", TU, &TrkUtil::funcNcl, 0.5, 1000., 1);
	f_ncl->SetNpx(500);
	f_ncl->Draw();
}
