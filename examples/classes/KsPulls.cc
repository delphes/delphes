#include "KsPulls.h"

//
// Final state configuration description
//
// Constructor
KsPulls::KsPulls()
{
	//
	// Pulls
	//
	// Pions
	// Parameters
	h_Pi_Par.push_back( new TH1D("h_Pi_D", "Pi impact parameter pull", 100, -10., 10.));
	h_Pi_Par.push_back( new TH1D("h_Pi_Phi0", "Pi Phi0 pull", 100, -10., 10.));
	h_Pi_Par.push_back( new TH1D("h_Pi_C", "Pi C pull", 100, -10., 10.));
	h_Pi_Par.push_back( new TH1D("h_Pi_Z0", "Pi Z0 pull", 100, -10., 10.));
	h_Pi_Par.push_back( new TH1D("h_Pi_Lm", "Pi Lm pull", 100, -10., 10.));
	// Momentum
	h_PiPx = new TH1D("h_PiPx","Pi X-momentum pull",100,-10.,10.);
	h_PiPy = new TH1D("h_PiPy","Pi Y-momentum pull",100,-10.,10.);
	h_PiPz = new TH1D("h_PiPz","Pi Z-momentum pull",100,-10.,10.);
	// Correlations
	h_PullPxCorrR = new TH2D("h_PullPxCorrR","Px pull vs R", 100, 0., 2000., 100, 0., 15.);
	h_PullPyCorrR = new TH2D("h_PullPyCorrR","Py pull vs R", 100, 0., 2000., 100, 0., 15.);
	h_PullPxCorrPt = new TH2D("h_PullPxCorrPt","Px pull vs Pt", 100, 0., 10., 100, 0., 15.);
	h_PullPyCorrPt = new TH2D("h_PullPyCorrPt","Py pull vs Pt", 100, 0., 10., 100, 0., 15.);
	h_PullPxCorrLm = new TH2D("h_PullPxCorrLm","Px pull vs {#lambda}", 100, -2.5, 2.5, 100, 0., 15.);
	h_PullPyCorrLm = new TH2D("h_PullPyCorrLm","Py pull vs {#lambda}", 100, -2.5, 2.5, 100, 0., 15.);
	h_PullPxCorrRD = new TH2D("h_PullPxCorrRD","Px pull vs R{^2}-D{^2}", 1000, -100., 400., 30, 0., 15.);
	h_PullPyCorrRD = new TH2D("h_PullPyCorrRD","Py pull vs R{^2}-D{^2}", 1000, -100., 400., 30, 0., 15.);
	// Ks vertex
	h_KsXv = new TH1D("h_KsXv","Ks X-vertex pull",100,-10.,10.);
	h_KsYv = new TH1D("h_KsYv","Ks Y-vertex pull",100,-10.,10.);
	h_KsZv = new TH1D("h_KsZv","Ks Z-vertex pull",100,-10.,10.);
	// Ks vertex distributions
	h_KsRgen = new TH1D("h_KsRgen","Generated Ks decay radius", 100., 0.,2000.);
	h_KsRrec = new TH1D("h_KsRrec","Reconstructe Ks decay radius", 100., 0.,2000.);
	TH1D* h_KsRrec;
	// Ks vertex momentum
	h_KsPx = new TH1D("h_KsPx","Ks X-momentum pull",100,-10.,10.);
	h_KsPy = new TH1D("h_KsPy","Ks Y-momentum pull",100,-10.,10.);
	h_KsPz = new TH1D("h_KsPz","Ks Z-momentum pull",100,-10.,10.);
	//
	// Ks masses
	h_KsdGenMass = new TH1D("h_KsdGenMass", "Ks generated mass difference", 100, -.05, .05);
	h_KsdRecMass = new TH1D("h_KsdRecMass", "Ks reconstructed mass difference", 100, -.05, .05);
	h_KsMassErr  = new TH1D("h_KsMassErr" , "Ks reconstructed mass error", 100, 0., 0.08);
	h_KsMassPull = new TH1D("h_KsMassPull", "Ks mass pull", 100, -10., 10.);
	//
	// Ks parameters
	h_Ks_Par.push_back( new TH1D("h_Ks_D", "Ks impact parameter pull", 100, -10., 10.));
	h_Ks_Par.push_back( new TH1D("h_Ks_Phi0", "Ks Phi0 pull", 100, -10., 10.));
	h_Ks_Par.push_back( new TH1D("h_Ks_Pt", "Ks Pt pull", 100, -10., 10.));
	h_Ks_Par.push_back( new TH1D("h_Ks_Z0", "Ks Z0 pull", 100, -10., 10.));
	h_Ks_Par.push_back( new TH1D("h_Ks_Lm", "Ks Lm pull", 100, -10., 10.));
}

//
// Destructor
KsPulls::~KsPulls() 
{
}

void KsPulls::Fill(VState* KsState, VertexMore* vKs)
{
	const Double_t Bz = 2.;
	const Int_t nPi = 2;
	GenParticle* pKs   = KsState->GetMother();	// Pointer to Ks genparticle
	//
	// Pions after fit
	//
	VertexFit* vvKs = vKs->GetVinput();		// Get pointer to vertex fit
	for(Int_t i=0; i<nPi; i++){			// Loop on pions in vertex
		GenParticle* Pion = KsState->GetGen(i);	// Pi one from  Ks
		TVector3 gvPi(Pion->X,  Pion->Y,  Pion->Z);	// Origin vertex
		TVectorD rvPi = vKs->GetXv();			// Ks reconstructed vertex
		Double_t R2 = rvPi(0)*rvPi(0)+rvPi(1)*rvPi(1);
		h_KsRgen->Fill(gvPi.Pt());
		h_KsRrec->Fill(TMath::Sqrt(R2));
		// Momentum
		TVector3 rpPi = vKs->GetMomentum(i);		// Reconstructed momentum
		TVector3 gpPi(Pion->Px, Pion->Py, Pion->Pz);	// Generated momentum
		TMatrixDSym cpPi = vKs->GetMomentumC(i);	// Covariance matrix
		// Parameters
		Double_t Q = vKs->GetCharge(i);
		gvPi *= 1.0e-3;							// Change to meters
		TVectorD genPar = TrkUtil::XPtoPar(gvPi, gpPi, Q, Bz);
		gvPi *= 1.0e3;
		TVectorD gParPi = TrkUtil::ParToMm(genPar);			// Back to mm
		TVectorD rParPi = vvKs->GetNewPar(i);				// Reconstructed
		Double_t D2 = rParPi(0)*rParPi(0);
		TMatrixDSym parPiCov = vvKs->GetNewCov(i);
		for(Int_t j=0; j<5; j++) h_Pi_Par[j]->Fill((rParPi(j)-gParPi(j))/TMath::Sqrt(parPiCov(j,j)));
		//
		Double_t PullPx = (rpPi(0)-gpPi(0))/TMath::Sqrt(cpPi(0,0));
		Double_t PullPy = (rpPi(1)-gpPi(1))/TMath::Sqrt(cpPi(1,1));
		Double_t PullPz = (rpPi(2)-gpPi(2))/TMath::Sqrt(cpPi(2,2));
		h_PiPx->Fill(PullPx);
		h_PiPy->Fill(PullPy);
		h_PiPz->Fill(PullPz);
		// Correlations
		Double_t limit = 4.0;
		if(TMath::Abs(PullPx)> limit){
			h_PullPxCorrR->Fill(gvPi.Pt(), TMath::Abs(PullPx));
			h_PullPyCorrR->Fill(gvPi.Pt(), TMath::Abs(PullPy));
			h_PullPxCorrPt->Fill(gpPi.Pt(), TMath::Abs(PullPx));
			h_PullPyCorrPt->Fill(gpPi.Pt(), TMath::Abs(PullPy));
			h_PullPxCorrRD->Fill(R2-D2, TMath::Abs(PullPx));
			h_PullPyCorrRD->Fill(R2-D2, TMath::Abs(PullPy));
			h_PullPxCorrLm->Fill(gParPi(4), TMath::Abs(PullPx));
			h_PullPyCorrLm->Fill(gParPi(4), TMath::Abs(PullPy));
		}
	}
	//
	//
	// Ks Vertex
	GenParticle* Pion1  = KsState->GetGen(0);	// Pi one from  Ks
	GenParticle* Pion2  = KsState->GetGen(1);	// Pi two from  Ks
	TVector3 gvKs(Pion1->X, Pion1->Y, Pion1->Z);	// K0s generated vertex
	TVectorD rvKs = vKs->GetXv();			// Ks reconstructed vertex
	TMatrixDSym cvKs = vKs->GetXvCov();		// Ks vertex covariance
	// Fill
	h_KsXv->Fill((rvKs(0)-gvKs(0))/TMath::Sqrt(cvKs(0,0)));
	h_KsYv->Fill((rvKs(1)-gvKs(1))/TMath::Sqrt(cvKs(1,1)));
	h_KsZv->Fill((rvKs(2)-gvKs(2))/TMath::Sqrt(cvKs(2,2)));
	// Ks Momentum
	TVector3 gpKs(pKs->Px, pKs->Py, pKs->Pz);	// Ks generated momentum
	TVector3 rpKs = vKs->GetTotalP();		// Ks reconstructed momentum
	TMatrixDSym cpKs = vKs->GetTotalPcov();		// Ks momentum covariance
	// Fill
	h_KsPx->Fill((rpKs(0)-gpKs(0))/TMath::Sqrt(cpKs(0,0)));
	h_KsPy->Fill((rpKs(1)-gpKs(1))/TMath::Sqrt(cpKs(1,1)));
	h_KsPz->Fill((rpKs(2)-gpKs(2))/TMath::Sqrt(cpKs(2,2)));
	//
	// Masses
	// Generated
	Double_t KsMass = pKs->Mass;
	Double_t PiMass = Pion1->Mass;
	TVector3 gpPi1(Pion1->Px, Pion1->Py, Pion1->Pz);
	Double_t gEPi1 = TMath::Sqrt(PiMass*PiMass+gpPi1.Mag2());
	TLorentzVector glPi1(gpPi1, gEPi1);
	//	
	TVector3 gpPi2(Pion2->Px, Pion2->Py, Pion2->Pz);
	Double_t gEPi2 = TMath::Sqrt(PiMass*PiMass+gpPi2.Mag2());
	TLorentzVector glPi2(gpPi2, gEPi2);
	TLorentzVector KsMomG = glPi1+glPi2;
	//
	h_KsdGenMass->Fill(KsMomG.M()-KsMass);
	//
	// Reconstructed
	TVector3 rpPi1 = vKs->GetMomentum(0);
	Double_t rEPi1 = TMath::Sqrt(PiMass*PiMass+rpPi1.Mag2());
	TLorentzVector rlPi1(rpPi1, rEPi1);
	//	
	TVector3 rpPi2 = vKs->GetMomentum(1);;
	Double_t rEPi2 = TMath::Sqrt(PiMass*PiMass+rpPi2.Mag2());
	TLorentzVector rlPi2(rpPi2, rEPi2);
	TLorentzVector KsMomR = rlPi1+rlPi2;
	//
	h_KsdRecMass->Fill(KsMomR.M()-KsMass);
	// Mass error
	Double_t Etot = KsMomG.E();
	TVector3 Ptot = KsMomG.Vect();
	TVector3 dm1 = (1./KsMass)*((Etot/gEPi1)*gpPi1-Ptot);	// dm/dp1
	TVector3 dm2 = (1./KsMass)*((Etot/gEPi2)*gpPi2-Ptot);	// dm/dp2
	Double_t dmp[6] = {dm1(0), dm1(1), dm1(2),dm2(0), dm2(1), dm2(2)};
	TVectorD DmDp(6, dmp);
	TMatrixDSym Kscov = vKs->GetBigCov();			
	TMatrixDSym pcov = Kscov.GetSub(3,8,3,8);		// cov(p1,p2)
	Double_t mKsErr = TMath::Sqrt(pcov.Similarity(DmDp));	// Error on Ks mass
	//
	h_KsMassErr->Fill(mKsErr);
	Double_t PullM = (KsMomR.M()-KsMass)/mKsErr;
	if(gvKs.Mag() < 100.) h_KsMassPull->Fill(PullM);
	if(TMath::Abs(PullM) > 4.0){
		std::cout<<"Large mass pull!"<<std::endl;
		std::cout<<"x rec: "<<rvKs(0)<<", "<<rvKs(1)<<", "<<rvKs(2)<<
			   ", p rec: "<<rpKs(0)<<", "<<rpKs(1)<<", "<<rpKs(2)<<std::endl;
		std::cout<<"x gen: "<<gvKs(0)<<", "<<gvKs(1)<<", "<<gvKs(2)<<
			   ", p gen: "<<gpKs(0)<<", "<<gpKs(1)<<", "<<gpKs(2)<<std::endl;
	}
	// 
	// Ks Parameters
	TVectorD gPar = TrkUtil::XPtoPar_N(gvKs, gpKs);		// Generated
	TVectorD rPar = vKs->GetVpar();				// Reconstructed
	TMatrixDSym parCov = vKs->GetVcov();
	for(Int_t i=0; i<5; i++) h_Ks_Par[i]->Fill((rPar(i)-gPar(i))/TMath::Sqrt(parCov(i,i)));
}

// Print 
void KsPulls::Print() 
{

	//
	// Ks plots
	TCanvas * cnv10 = new TCanvas("cnv10","Ks vertex pulls",10,10, 900,600);
	cnv10->Divide(3,2);
	// X
	cnv10->cd(1); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_KsXv->Fit("gaus"); h_KsXv->Draw();
	// Y
	cnv10->cd(2); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_KsYv->Fit("gaus"); h_KsYv->Draw();
	// Z
	cnv10->cd(3); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_KsZv->Fit("gaus"); h_KsZv->Draw();
	// Px
	cnv10->cd(4); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_KsPx->Fit("gaus"); h_KsPx->Draw();
	// Py
	cnv10->cd(5); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_KsPy->Fit("gaus"); h_KsPy->Draw();
	// Pz
	cnv10->cd(6); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_KsPz->Fit("gaus"); h_KsPz->Draw();
	//
	// Mass
	TCanvas * cnv11 = new TCanvas("cnv11","Ks generated mass difference",100,100, 900,600);
	cnv11->Divide(2,2);
	cnv11->cd(1);
	h_KsdGenMass->Draw();
	cnv11->cd(2);
	h_KsdRecMass->Draw();
	cnv11->cd(3);
	h_KsMassErr->Draw();
	cnv11->cd(4); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_KsMassPull->Fit("gaus"); h_KsMassPull->Draw();
	//
	// Parameters
	TCanvas * cnv12 = new TCanvas("cnv12","Ks track parameter pulls",100,100, 900,600);
	cnv12->Divide(3,2);
	for(Int_t i=0; i<5; i++){
		cnv12->cd(i+1); gPad->SetLogy(1);
		gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
		h_Ks_Par[i]->Fit("gaus"); h_Ks_Par[i]->Draw();
	}
	//
	// Pion plots
	// Parameters
	TCanvas * cnv13 = new TCanvas("cnv13","Pi track parameter pulls",150,150, 900,600);
	cnv13->Divide(3,2);
	for(Int_t i=0; i<5; i++){
		cnv13->cd(i+1); gPad->SetLogy(1);
		gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
		h_Pi_Par[i]->Fit("gaus"); h_Pi_Par[i]->Draw();
	}
	// Momentum
	// Px
	TCanvas * cnv14 = new TCanvas("cnv14","Pi momentum pulls",200,200, 900,600);
	cnv14->Divide(3,1);
	cnv14->cd(1); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_PiPx->Fit("gaus"); h_PiPx->Draw();
	// Py
	cnv14->cd(2); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_PiPy->Fit("gaus"); h_PiPy->Draw();
	// Pz
	cnv14->cd(3); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_PiPz->Fit("gaus"); h_PiPz->Draw();
	// Correlations
	TCanvas * cnv15 = new TCanvas("cnv15","Pi momentum pull correlations ",250,250, 900,600);
	cnv15->Divide(2,2);
	cnv15->cd(1); h_PullPxCorrR-> SetMarkerStyle(3); h_PullPxCorrR->Draw();
	cnv15->cd(2); 
	h_KsRgen->Draw(); h_KsRrec->SetLineColor(kRed); h_KsRrec->Draw("SAME"); 
	cnv15->cd(3); h_PullPxCorrPt-> SetMarkerStyle(3); h_PullPxCorrPt->Draw();
	cnv15->cd(4); h_PullPxCorrLm-> SetMarkerStyle(3); h_PullPxCorrLm->Draw();
}
