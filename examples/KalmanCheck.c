
#include <iostream>
#include <TString.h>
#include <TLegend.h>
#include "examples/classes/KalmanCk.h"
//
//
void KalmanCheck(TString GEOM, Bool_t Old=kFALSE, Bool_t Res=kTRUE, Bool_t MS=kTRUE)
{
//
// Initialize geometry
	std::cout << "Text geometry file: " << GEOM << std::endl;
	char* GeoName = (char*) GEOM.Data();
	SolGeom* G = new SolGeom(GeoName);
	// Start plotting class
	KalmanCk * K = new KalmanCk(G);
	K->DrawPtScan(3);	// Draw sample tracks
	K->SetMode(Res, MS);	// Set tracking mode
	K->SetOld(Old);	// Select Standard calculation version
	// Fill graphs
	K->Fill();
	// Display
	K->Print();
	//
	// Publication plots
	//
	/*
	TCanvas *cs = new TCanvas("cs","Standard method - D resolutions",200,200, 900, 500);
	cs->Divide(1,1);
	TCanvas *ck = new TCanvas("ck","Kalman method - D resolutions",250,250, 900, 500);
	ck->Divide(1,1);
	//
	TLegend* lg = new TLegend(0.2, 0.9, 0.6, 0.70); 
	TString *Langa = new TString[K->fNptFpt];
	// loop over ref pt
	for(Int_t np=0; np<K->fNptFpt; np++){
		//
		
		Langa[np].Form("Track momentum %.0f GeV", K->fPtFpt[np]);
		lg->AddEntry(K->gs_D_Ang[np], Langa[np], "L");
		//
		TString axis = "#theta (degrees)";
		// D plots vs Angle
		cs->cd(1); gPad->SetLogy(1);	// standard D plots
		TString title_D = "#sigma(D) #mum";
		Double_t ymax_D = 100.;
		Double_t ymin_D = 1.0;
		K->GrSetup(K->gs_D_Ang[np],ymax_D,np+1,title_D,axis);
		K->gs_D_Ang[np]->SetMinimum(1.0);
		lg->Draw();
		if(np == 0)K->gs_D_Ang[np]->Draw("APL");
		else K->gs_D_Ang[np]->Draw("SAMEPL");
		ck->cd(1); gPad->SetLogy(1);	// Kalman D plots
		K->GrSetup(K->gk_D_Ang[np],ymax_D,np+1,title_D,axis);
		K->gk_D_Ang[np]->SetMinimum(1.0);
		lg->Draw();
		if(np == 0)K->gk_D_Ang[np]->Draw("APL");
		else K->gk_D_Ang[np]->Draw("SAMEPL");
	}
	*/
}

