#include "KalmanCk.h"

//
// Compare Kalman and full fit tracking errors
//
// Constructor
KalmanCk::KalmanCk(SolGeom *G)
{
//
// Initialize geometry
//
	fG = G;				// Initialize geometry
	fBz = fG->B();			// Magnetic field
	//
	// Operation modes
	fRes = kTRUE;	// Turn on detector resolution
	fMS  = kTRUE;	// Turn on multiple scattering
	fOld = kFALSE;  // Assume new method is used
// Pt scans at fixed angles
//
	fPt_fixA =  {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
		     1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 
		     10., 20., 30., 40., 50., 60., 70., 80., 90.,
		     100.,200.};
	fNpt_Fa = (Int_t)fPt_fixA.size();
	
	fAngFa = {10., 20., 30., 60., 90.};
	fNangFa = (Int_t)fAngFa.size();
//
// Polar angle  scans at fixed pt
//
	fAng_fixPt = {10., 15., 20., 25., 30., 35., 40.,
		      45., 50., 55., 60., 65., 70., 75.,
		      80., 85., 90.};
	fNang_Fpt = (Int_t)fAng_fixPt.size();
	//
	fPtFpt = {1., 10., 100.};
	fNptFpt = (Int_t)fPtFpt.size();
}

//
// Destructor
KalmanCk::~KalmanCk() 
{
	// Clear all configuration vectors
	fPt_fixA.clear();
	fAngFa.clear();
	fAng_fixPt.clear();
	fPtFpt.clear();
	// Clear all graphs
	// Standard resolution graphs for pt scans
	gs_D_Pt.clear();	// D resolution graphs for pt scans
	gs_Phi0_Pt.clear();// Phi0 resolution graphs for pt scans
	gs_Pt_Pt.clear();	// Pt resolution graphs for pt scans
	gs_z0_Pt.clear();	// z0 resolution graphs for pt scans
	gs_cot_Pt.clear();	// cot(theta) graphs for pt scans
	//
	// Kalman resolution graphs for pt scans
	gk_D_Pt.clear();	// D resolution graphs for pt scans
	gk_Phi0_Pt.clear();// Phi0 resolution graphs for pt scans
	gk_Pt_Pt.clear();	// Pt resolution graphs for pt scans
	gk_z0_Pt.clear();	// z0 resolution graphs for pt scans
	gk_cot_Pt.clear();	// cot(theta) graphs for pt scans
	//
	// Standard resolution graphs for angle scans
	gs_D_Ang.clear();		// D resolution graphs for angle scans
	gs_Phi0_Ang.clear();	// Phi0 resolution graphs for angle scans
	gs_Pt_Ang.clear();		// Pt resolution graphs for angle scans
	gs_z0_Ang.clear();		// z0 resolution graphs for angle scans
	gs_cot_Ang.clear();	// cot(theta) graphs for angle scans
	//
	// Kalman resolution graphs for angle scans
	gk_D_Ang.clear();		// D resolution graphs for angle scans
	gk_Phi0_Ang.clear();	// Phi0 resolution graphs for angle scans
	gk_Pt_Ang.clear();		// Pt resolution graphs for angle scans
	gk_z0_Ang.clear();		// z0 resolution graphs for angle scans
	gk_cot_Ang.clear();	// cot(theta) graphs for angle scans
}

void KalmanCk::SetMode(Bool_t Res, Bool_t MS)
{
	fRes = Res;
	fMS  = MS;
	if(fRes == kFALSE && fMS == kFALSE){
		std::cout<<
		"Setting both Res and MS to FALSE is not allowed. Reset both to TRUE"
		<<std::endl;
		fRes = kTRUE;
		fMS  = kTRUE;
	}	
}
void KalmanCk::Fill()
{
//==============================================================================
//
// Fill Pt scan plots
//
//==============================================================================
//
//	Standard
	std::vector<TVectorD*>vsD_Pt;
	std::vector<TVectorD*>vsPhi0_Pt;
	std::vector<TVectorD*>vsPt_Pt;
	std::vector<TVectorD*>vsZ0_Pt;
	std::vector<TVectorD*>vsCot_Pt;
//	Kalman
	std::vector<TVectorD*>vkD_Pt;
	std::vector<TVectorD*>vkPhi0_Pt;
	std::vector<TVectorD*>vkPt_Pt;
	std::vector<TVectorD*>vkZ0_Pt;
	std::vector<TVectorD*>vkCot_Pt;
	//
	// Start scanning angles and transverse momenta
	//
	for(Int_t na=0; na<fNangFa; na++){	// Angle loop
		Double_t angle = fAngFa[na];	// Track angle in degrees
		Double_t theta = angle * TMath::Pi() / 180.;	// Convert to radians
		// Standard resolution arrays
		TVectorD asD_Pt   (fNpt_Fa);	// Transverse impact parameter
		TVectorD asPhi0_Pt(fNpt_Fa);	// Phi direction at min approach
		TVectorD asPt_Pt  (fNpt_Fa);	// transverse momentum
		TVectorD asZ0_Pt  (fNpt_Fa);	// Z at min. approach
		TVectorD asCot_Pt (fNpt_Fa);	// Cotangent of polar angle
		// Kalman resolution arrays
		TVectorD akD_Pt   (fNpt_Fa);
		TVectorD akPhi0_Pt(fNpt_Fa);
		TVectorD akPt_Pt  (fNpt_Fa);
		TVectorD akZ0_Pt  (fNpt_Fa);
		TVectorD akCot_Pt (fNpt_Fa);
		//
		for(Int_t np=0; np<fNpt_Fa; np++){	// Pt loop
			Double_t pt = fPt_fixA[np];	// track pt
			Double_t pz = 0;
			if (angle != 90.) pz = pt / TMath::Tan(theta);
			//
			TVector3 p(pt, 0.0, pz);	// Track 3-momentum
			TVector3 x(0.0, 0.0, 0.0);	// Track starting point
			// Standard results
			SolTrack Strack(x, p, fG);	// Initialize track
			if(fOld)Strack.OldCovCalc(fRes,fMS);
			else Strack.CovCalc(fRes,fMS);
			asD_Pt(np)    = Strack.s_D()*1.e6;	// Convert to microns
			asPhi0_Pt(np) = Strack.s_phi0();
			asPt_Pt(np)   = Strack.s_pt();		// sigma(pt)/pt
			asZ0_Pt(np)   = Strack.s_z0()*1.e6;	// Convert to microns
			asCot_Pt(np)  = Strack.s_ct();
			// Kalman results
			SolTrack Ktrack(x, p, fG);	// Initialize track
			Ktrack.KalmanCov(fRes, fMS);
			akD_Pt(np)    = Ktrack.s_D()*1.e6;	// Convert to microns
			akPhi0_Pt(np) = Ktrack.s_phi0();
			akPt_Pt(np)   = Ktrack.s_pt();		// sigma(pt)/pt
			akZ0_Pt(np)   = Ktrack.s_z0()*1.e6;	// Convert to microns
			akCot_Pt(np)  = Ktrack.s_ct();
		}
		//
		// Store results and fill graphs
		//
		// Standard arrays
		vsD_Pt.push_back(new TVectorD(asD_Pt));
		vsPhi0_Pt.push_back(new TVectorD(asPhi0_Pt));
		vsPt_Pt.push_back(	new TVectorD(asPt_Pt));
		vsZ0_Pt.push_back(	new TVectorD(asZ0_Pt));
		vsCot_Pt.push_back(	new TVectorD(asCot_Pt));
		//
		// Standard graphs
		gs_D_Pt.push_back(new TGraph(fNpt_Fa,fPt_fixA.data(),
							vsD_Pt[na]->GetMatrixArray()));
		gs_Phi0_Pt.push_back(new TGraph(fNpt_Fa,fPt_fixA.data(),
					        vsPhi0_Pt[na]->GetMatrixArray()));
		gs_Pt_Pt.push_back(new TGraph(fNpt_Fa,fPt_fixA.data(),
							vsPt_Pt[na]->GetMatrixArray()));
		gs_z0_Pt.push_back(new TGraph(fNpt_Fa,fPt_fixA.data(),
							vsZ0_Pt[na]->GetMatrixArray()));
		gs_cot_Pt.push_back(new TGraph(fNpt_Fa,fPt_fixA.data(),
							vsCot_Pt[na]->GetMatrixArray()));
		//
		// Kalman arrays
		vkD_Pt.push_back(	new TVectorD(akD_Pt));
		vkPhi0_Pt.push_back(new TVectorD(akPhi0_Pt));
		vkPt_Pt.push_back(	new TVectorD(akPt_Pt));
		vkZ0_Pt.push_back(	new TVectorD(akZ0_Pt));
		vkCot_Pt.push_back(	new TVectorD(akCot_Pt));
		//
		// Kalman graphs
		gk_D_Pt.push_back(new TGraph(fNpt_Fa,fPt_fixA.data(),
							vkD_Pt[na]->GetMatrixArray()));
		gk_Phi0_Pt.push_back(new TGraph(fNpt_Fa,fPt_fixA.data(),
					        vkPhi0_Pt[na]->GetMatrixArray()));
		gk_Pt_Pt.push_back(new TGraph(fNpt_Fa,fPt_fixA.data(),
							vkPt_Pt[na]->GetMatrixArray()));
		gk_z0_Pt.push_back(new TGraph(fNpt_Fa,fPt_fixA.data(),
							vkZ0_Pt[na]->GetMatrixArray()));
		gk_cot_Pt.push_back(new TGraph(fNpt_Fa,fPt_fixA.data(),
							vkCot_Pt[na]->GetMatrixArray()));
	}
//
//==============================================================================
//
// Fill angle scan plots
//
//==============================================================================
//
//	Standard
	std::vector<TVectorD*>vsD_Ang;
	std::vector<TVectorD*>vsPhi0_Ang;
	std::vector<TVectorD*>vsPt_Ang;
	std::vector<TVectorD*>vsZ0_Ang;
	std::vector<TVectorD*>vsCot_Ang;
//	Kalman
	std::vector<TVectorD*>vkD_Ang;
	std::vector<TVectorD*>vkPhi0_Ang;
	std::vector<TVectorD*>vkPt_Ang;
	std::vector<TVectorD*>vkZ0_Ang;
	std::vector<TVectorD*>vkCot_Ang;
	//
	// Start scanning angles and transverse momenta
	//
	for(Int_t np=0; np<fNptFpt; np++){	// Pt loop
		// Standard resolution arrays
		TVectorD asD_Ang   (fNang_Fpt);	// Transverse impact parameter
		TVectorD asPhi0_Ang(fNang_Fpt);	// Phi direction at min approach
		TVectorD asPt_Ang  (fNang_Fpt);	// transverse momentum
		TVectorD asZ0_Ang  (fNang_Fpt);	// Z at min. approach
		TVectorD asCot_Ang (fNang_Fpt);	// Cotangent of polar angle
		// Kalman resolution arrays
		TVectorD akD_Ang   (fNang_Fpt);
		TVectorD akPhi0_Ang(fNang_Fpt);
		TVectorD akPt_Ang  (fNang_Fpt);
		TVectorD akZ0_Ang  (fNang_Fpt);
		TVectorD akCot_Ang (fNang_Fpt);
		//
		for(Int_t na=0; na<fNang_Fpt; na++){	// Angle loop
			Double_t angle = fAng_fixPt[na];	// Track angle in degrees
			Double_t theta = angle * TMath::Pi() / 180.;	// Convert to radians
			Double_t pt = fPtFpt[np];	// track pt
			Double_t pz = 0;
			if (angle != 90.) pz = pt / TMath::Tan(theta);
			//
			TVector3 p(pt, 0.0, pz);	// Track 3-momentum
			TVector3 x(0.0, 0.0, 0.0);	// Track starting point
			// Standard results
			SolTrack Strack(x, p, fG);	// Initialize track
			if(fOld)Strack.OldCovCalc(fRes,fMS);
			else Strack.CovCalc(fRes,fMS);
			asD_Ang(na)    = Strack.s_D()*1.e6;	// Convert to microns
			asPhi0_Ang(na) = Strack.s_phi0();
			asPt_Ang(na)   = Strack.s_pt();		// sigma(pt)/pt
			asZ0_Ang(na)   = Strack.s_z0()*1.e6;	// Convert to microns
			asCot_Ang(na)  = Strack.s_ct();
			// Kalman results
			SolTrack Ktrack(x, p, fG);	// Initialize track
			Ktrack.KalmanCov(fRes, fMS);
			akD_Ang(na)    = Ktrack.s_D()*1.e6;	// Convert to microns
			akPhi0_Ang(na) = Ktrack.s_phi0();
			akPt_Ang(na)   = Ktrack.s_pt();		// sigma(pt)/pt
			akZ0_Ang(na)   = Ktrack.s_z0()*1.e6;	// Convert to microns
			akCot_Ang(na)  = Ktrack.s_ct();
		}
		//
		// Store results and fill graphs
		//
		// Standard arrays
		vsD_Ang.push_back(new TVectorD(asD_Ang));
		vsPhi0_Ang.push_back(new TVectorD(asPhi0_Ang));
		vsPt_Ang.push_back(new TVectorD(asPt_Ang));
		vsZ0_Ang.push_back(new TVectorD(asZ0_Ang));
		vsCot_Ang.push_back(new TVectorD(asCot_Ang));
		//
		// Standard graphs
		gs_D_Ang.push_back(new TGraph(fNang_Fpt,fAng_fixPt.data(),
		                             vsD_Ang[np]->GetMatrixArray()));
		gs_Phi0_Ang.push_back(new TGraph(fNang_Fpt,fAng_fixPt.data(),
					        vsPhi0_Ang[np]->GetMatrixArray()));
		gs_Pt_Ang.push_back(new TGraph(fNang_Fpt,fAng_fixPt.data(),
				   	      vsPt_Ang[np]->GetMatrixArray()));
		gs_z0_Ang.push_back(new TGraph(fNang_Fpt,fAng_fixPt.data(),
					      vsZ0_Ang[np]->GetMatrixArray()));
		gs_cot_Ang.push_back(new TGraph(fNang_Fpt,fAng_fixPt.data(),
					       vsCot_Ang[np]->GetMatrixArray()));
		//
		// Kalman arrays
		vkD_Ang.push_back(new TVectorD(akD_Ang));
		vkPhi0_Ang.push_back(new TVectorD(akPhi0_Ang));
		vkPt_Ang.push_back(new TVectorD(akPt_Ang));
		vkZ0_Ang.push_back(new TVectorD(akZ0_Ang));
		vkCot_Ang.push_back(new TVectorD(akCot_Ang));
		//
		// Kalman graphs
		gk_D_Ang.push_back(new TGraph(fNang_Fpt,fAng_fixPt.data(),
		                             vkD_Ang[np]->GetMatrixArray()));
		gk_Phi0_Ang.push_back(new TGraph(fNang_Fpt,fAng_fixPt.data(),
					        vkPhi0_Ang[np]->GetMatrixArray()));
		gk_Pt_Ang.push_back(new TGraph(fNang_Fpt,fAng_fixPt.data(),
				   	      vkPt_Ang[np]->GetMatrixArray()));
		gk_z0_Ang.push_back(new TGraph(fNang_Fpt,fAng_fixPt.data(),
					      vkZ0_Ang[np]->GetMatrixArray()));
		gk_cot_Ang.push_back(new TGraph(fNang_Fpt,fAng_fixPt.data(),
					       vkCot_Ang[np]->GetMatrixArray()));
	}
}
//
// Display
void KalmanCk::DrawPtScan(Int_t i){
	//
	// Display canvas
	fG->Draw();
	TCanvas *Cnv = fG->cnv();
	Cnv->Draw();
	//
	// Loop over angles
	//
	if(i<0 || i>= fNpt_Fa){
		std::cout<<"KalmanCh::DrawPtScan: Index out of range. Set to "
		<<fNpt_Fa-1<<std::endl;}
		
	Double_t pt = fPt_fixA[i];
	TGraph **gr_trk = new TGraph*[fNangFa];
	for(Int_t i=0; i<fNangFa; i++){
		// Initialize track
		Double_t angle = fAngFa[i];
		Double_t theta = angle * TMath::Pi() / 180.;
		Double_t pz = 0;
		if (angle != 90.) pz = pt / TMath::Tan(theta);
		TVector3 p(pt, 0.0, pz);
		TVector3 x(0.0, 0.0, 0.0);
		SolTrack* trk_id = new SolTrack(x, p, fG);
		gr_trk[i] = trk_id->TrkPlot();	// graph intersection with layers
		gr_trk[i]->Draw("PLSAME");
	}
}
// Print 

void KalmanCk::GrSetup(TGraph* &g, Double_t yMax, Int_t Color, TString title, TString axis) 
{
	g->SetLineColor(Color);
	g->SetMarkerColor(Color);
	g->SetTitle(title);
	g->SetMinimum(0.0);
	if (yMax != 0.0) g->SetMaximum(yMax);
	g->GetXaxis()->SetTitle(axis);
}
void KalmanCk::Print() 
{
//==============================================================================
//
// Display momentum scans
//
//==============================================================================
//
	TCanvas *csPt;
	if(fOld){
		csPt = new TCanvas("csPt","OLD Standard method - Pt scans",50,50, 900, 500);
	}
	else {
		csPt = new TCanvas("csPt","Standard method - Pt scans",50,50, 900, 500);
	}
	csPt->Divide(3,2);
	TCanvas *ckPt = new TCanvas("ckPt","Kalman method - Pt scans",150,150, 900, 500);
	ckPt->Divide(3,2);
	//
	// Legends
	//
	TLegend* lsPt = new TLegend(0.2, 0.9, 0.6, 0.70);
	TString *LsPt = new TString[fNangFa];
	//
	TLegend* lkPt = new TLegend(0.2, 0.9, 0.6, 0.70);
	TString *LkPt = new TString[fNangFa];
	//
	for(Int_t na=0; na<fNangFa; na++){
		//
		LsPt[na].Form("Track angle %.0f#circ", fAngFa[na]);
		LkPt[na].Form("Track angle %.0f#circ", fAngFa[na]);
		lsPt->AddEntry(gs_D_Pt[na], LsPt[na], "L");
		lkPt->AddEntry(gk_D_Pt[na], LkPt[na], "L");
	}
	//
	// loop over ref angles
	Int_t iColor = 0;
	for(Int_t na=0; na<fNangFa; na++){
		//
		TString axis = "p_{t} (GeV/c)";
		iColor++;
		if(iColor == 5 or iColor == 10)iColor++;	// Avoid yellow and white
		// D plots vs Pt
		csPt->cd(1); gPad->SetLogy(1);	// standard D plots
		TString title_D = "#sigma(D) #mum";
		Double_t ymax_D = 100.;
		GrSetup(gs_D_Pt[na],ymax_D,iColor,title_D,axis);
		if(na == 0)gs_D_Pt[na]->Draw("APL");
		else gs_D_Pt[na]->Draw("SAMEPL");
		lsPt->Draw("SAME");
		ckPt->cd(1); gPad->SetLogy(1);	// Kalman D plots
		GrSetup(gk_D_Pt[na],ymax_D,iColor,title_D,axis);
		if(na == 0)gk_D_Pt[na]->Draw("APL");
		else gk_D_Pt[na]->Draw("SAMEPL");
		lkPt->Draw("SAME");
		//
		// Phi0 plots vs Pt
		csPt->cd(2); gPad->SetLogy(1);	// standard Phi0 plots
		TString title_phi = "#sigma(#varphi_{0}) rad";
		Double_t ymax_phi = 25.e-3;
		GrSetup(gs_Phi0_Pt[na],ymax_phi,iColor,title_phi,axis);
		if(na == 0)gs_Phi0_Pt[na]->Draw("APL");
		else gs_Phi0_Pt[na]->Draw("SAMEPL");
		lsPt->Draw("SAME");
		ckPt->cd(2); gPad->SetLogy(1);	// Kalman Phi0 plots
		GrSetup(gk_Phi0_Pt[na],ymax_phi,iColor,title_phi,axis);
		if(na == 0)gk_Phi0_Pt[na]->Draw("APL");
		else gk_Phi0_Pt[na]->Draw("SAMEPL");
		lkPt->Draw("SAME");
		//
		// Pt plots vs Pt
		csPt->cd(3); gPad->SetLogy(1);	// standard Pt plots
		TString title_pt = "#sigma(p_{t})/p_{t}";
		Double_t ymax_pt = 4.0e-2;
		GrSetup(gs_Pt_Pt[na],ymax_pt,iColor,title_pt,axis);
		if(na == 0)gs_Pt_Pt[na]->Draw("APL");
		else gs_Pt_Pt[na]->Draw("SAMEPL");
		lsPt->Draw("SAME");
		ckPt->cd(3); gPad->SetLogy(1);	// Kalman Pt plots
		GrSetup(gk_Pt_Pt[na],ymax_pt,iColor,title_pt,axis);
		if(na == 0)gk_Pt_Pt[na]->Draw("APL");
		else gk_Pt_Pt[na]->Draw("SAMEPL");
		lkPt->Draw("SAME");
		//
		// Z0 plots vs Pt
		csPt->cd(4); gPad->SetLogy(1);	// standard Z0 plots
		TString title_z = "#sigma(z_{0}) #mum";
		Double_t ymax_z = 100.;
		GrSetup(gs_z0_Pt[na],ymax_z,iColor,title_z,axis);
		if(na == 0)gs_z0_Pt[na]->Draw("APL");
		else gs_z0_Pt[na]->Draw("SAMEPL");
		lsPt->Draw("SAME");
		ckPt->cd(4); gPad->SetLogy(1);	// Kalman Z0 plots
		GrSetup(gk_z0_Pt[na],ymax_z,iColor,title_z,axis);
		if(na == 0)gk_z0_Pt[na]->Draw("APL");
		else gk_z0_Pt[na]->Draw("SAMEPL");
		lkPt->Draw("SAME");
		//
		// Cot(theta) plots vs Pt
		csPt->cd(5); gPad->SetLogy(1);	// standard cot(theta) plots
		TString title_ct = "#sigma(ctg #theta)";
		Double_t ymax_ct = .2;
		gs_cot_Pt[na]->SetMinimum(1.0E-5);
		gs_cot_Pt[na]->GetYaxis()->SetLimits(1.0E-5,ymax_ct);
		GrSetup(gs_cot_Pt[na],ymax_ct,iColor,title_ct,axis);
		if(na == 0)gs_cot_Pt[na]->Draw("APL");
		else gs_cot_Pt[na]->Draw("SAMEPL");
		lsPt->Draw("SAME");
		ckPt->cd(5); gPad->SetLogy(1);	// Kalman cot(theta) plots
		gk_cot_Pt[na]->SetMinimum(1.0E-5);
		gk_cot_Pt[na]->GetYaxis()->SetLimits(1.0E-5,ymax_ct);
		GrSetup(gk_cot_Pt[na],ymax_ct,iColor,title_ct,axis);
		if(na == 0)gk_cot_Pt[na]->Draw("APL");
		else gk_cot_Pt[na]->Draw("SAMEPL");
		lkPt->Draw("SAME");
	}
//
//==============================================================================
//
// Display angle scans
//
//==============================================================================
//
	TCanvas *csAng;
	if(fOld){
		csAng = new TCanvas("csAng","OLD Standard method - Angle scans",200,200, 900, 500);
	}
	else{
		csAng = new TCanvas("csAng","Standard method - Angle scans",200,200, 900, 500);
	}
	csAng->Divide(3,2);
	TCanvas *ckAng = new TCanvas("ckAng","Kalman method - Angle scans",250,250, 900, 500);
	//
	// Legends
	//
	TLegend* lsAng = new TLegend(0.2, 0.9, 0.6, 0.70);
	TString *LsAng = new TString[fNptFpt];
	//
	TLegend* lkAng = new TLegend(0.2, 0.9, 0.6, 0.70);
	TString *LkAng = new TString[fNptFpt];
	//
	for(Int_t np=0; np<fNptFpt; np++){
		//
		LsAng[np].Form("Track pt %.1f GeV/c", fPtFpt[np]);
		LkAng[np].Form("Track pt %.1f GeV/c", fPtFpt[np]);
		lsAng->AddEntry(gs_D_Ang[np], LsAng[np], "L");
		lkAng->AddEntry(gk_D_Ang[np], LkAng[np], "L");
	}
	ckAng->Divide(3,2);
	//
	// loop over ref pt
	for(Int_t np=0; np<fNptFpt; np++){
		//
		TString axis = "#theta (degrees)";
		// D plots vs Angle
		csAng->cd(1); gPad->SetLogy(1);	// standard D plots
		TString title_D = "#sigma(D) #mum";
		Double_t ymax_D = 100.;
		GrSetup(gs_D_Ang[np],ymax_D,np+1,title_D,axis);
		if(np == 0)gs_D_Ang[np]->Draw("APL");
		else gs_D_Ang[np]->Draw("SAMEPL");
		lsAng->Draw("SAME");
		ckAng->cd(1); gPad->SetLogy(1);	// Kalman D plots
		GrSetup(gk_D_Ang[np],ymax_D,np+1,title_D,axis);
		if(np == 0)gk_D_Ang[np]->Draw("APL");
		else gk_D_Ang[np]->Draw("SAMEPL");
		lkAng->Draw("SAME");
		//
		// Phi0 plots vs Angle
		csAng->cd(2); gPad->SetLogy(1);	// standard Phi0 plots
		TString title_phi = "#sigma(#varphi_{0}) rad";
		Double_t ymax_phi = 25.e-3;
		GrSetup(gs_Phi0_Ang[np],ymax_phi,np+1,title_phi,axis);
		if(np == 0)gs_Phi0_Ang[np]->Draw("APL");
		else gs_Phi0_Ang[np]->Draw("SAMEPL");
		lsAng->Draw("SAME");
		ckAng->cd(2); gPad->SetLogy(1);	// Kalman Phi0 plots
		GrSetup(gk_Phi0_Ang[np],ymax_phi,np+1,title_phi,axis);
		if(np == 0)gk_Phi0_Ang[np]->Draw("APL");
		else gk_Phi0_Ang[np]->Draw("SAMEPL");
		lkAng->Draw("SAME");
		//
		// Pt plots vs Angle
		csAng->cd(3); gPad->SetLogy(1);	// standard Pt plots
		TString title_pt = "#sigma(p_{t})/p_{t}";
		Double_t ymax_pt = 4.0e-2;
		GrSetup(gs_Pt_Ang[np],ymax_pt,np+1,title_pt,axis);
		if(np == 0)gs_Pt_Ang[np]->Draw("APL");
		else gs_Pt_Ang[np]->Draw("SAMEPL");
		lsAng->Draw("SAME");
		ckAng->cd(3); gPad->SetLogy(1);	// Kalman Pt plots
		GrSetup(gk_Pt_Ang[np],ymax_pt,np+1,title_pt,axis);
		if(np == 0)gk_Pt_Ang[np]->Draw("APL");
		else gk_Pt_Ang[np]->Draw("SAMEPL");
		lkAng->Draw("SAME");
		//
		// Z0 plots vs Ang
		csAng->cd(4); gPad->SetLogy(1);	// standard Z0 plots
		TString title_z = "#sigma(z_{0}) #mum";
		Double_t ymax_z = 100.;
		GrSetup(gs_z0_Ang[np],ymax_z,np+1,title_z,axis);
		if(np == 0)gs_z0_Ang[np]->Draw("APL");
		else gs_z0_Ang[np]->Draw("SAMEPL");
		lsAng->Draw("SAME");
		ckAng->cd(4); gPad->SetLogy(1);	// Kalman Z0 plots
		GrSetup(gk_z0_Ang[np],ymax_z,np+1,title_z,axis);
		if(np == 0)gk_z0_Ang[np]->Draw("APL");
		else gk_z0_Ang[np]->Draw("SAMEPL");
		lkAng->Draw("SAME");
		//
		// Cot(theta) plots vs Angle
		csAng->cd(5); gPad->SetLogy(1);	// standard cot(theta) plots
		TString title_ct = "#sigma(ctg #theta)";
		Double_t ymax_ct = .2;
		gs_cot_Ang[np]->SetMinimum(1.0E-5);
		gs_cot_Ang[np]->GetYaxis()->SetLimits(1.0E-5,ymax_ct);
		GrSetup(gs_cot_Ang[np],ymax_ct,np+1,title_ct,axis);
		if(np == 0)gs_cot_Ang[np]->Draw("APL");
		else gs_cot_Ang[np]->Draw("SAMEPL");
		lsAng->Draw("SAME");
		ckAng->cd(5); gPad->SetLogy(1);	// Kalman cot(theta) plots
		gk_cot_Ang[np]->SetMinimum(1.0E-5);
		gk_cot_Ang[np]->GetYaxis()->SetLimits(1.0E-5,ymax_ct);
		GrSetup(gk_cot_Ang[np],ymax_ct,np+1,title_ct,axis);
		if(np == 0)gk_cot_Ang[np]->Draw("APL");
		else gk_cot_Ang[np]->Draw("SAMEPL");
		lkAng->Draw("SAME");
	}

	
}
