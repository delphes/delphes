#include "VtxPulls.h"

//
// Final state configuration description
//
// Constructor
VtxPulls::VtxPulls()
{
	//
	// Pulls
	// Ds vertex
	h_Xv = new TH1D("h_Xv","X-vertex pull",100,-10.,10.);
	h_Yv = new TH1D("h_Yv","Y-vertex pull",100,-10.,10.);
	h_Zv = new TH1D("h_Zv","Z-vertex pull",100,-10.,10.);
	h_Chi2 = new TH1D("h_Chi2", "Vertex #chi^{2}/N_{dof}", 100, 0., 10.);
	// Updated track parameter pulls after vertex fit
	h_D  = new TH1D("h_D","D impact parameter pull",100,-10.,10.);
	h_P0 = new TH1D("h_P0","Phi0 parameter pull",100,-10.,10.);
	h_C  = new TH1D("h_C","Half curvature pull",100,-10.,10.);
	h_Z0 = new TH1D("h_Z0","Z0 pull",100,-10.,10.);
	h_Ct = new TH1D("h_Ct","Cot theta pull",100,-10.,10.);
	// track momenta at vertex
	h_Px = new TH1D("h_Px","Momentum x pull",100,-10.,10.);
	h_Py = new TH1D("h_Py","Momentum y pull",100,-10.,10.);
	h_Pz = new TH1D("h_Pz","Momentum z pull",100,-10.,10.);
}

//
// Destructor
VtxPulls::~VtxPulls() 
{
}

void VtxPulls::Fill(TVectorD vTrue, std::vector<TVectorD> pTrue, VertexFit *Vfit)
{
	// Global variables
	Double_t Bz = 2.;	// Magnetic field
	//
	// Vertex pulls
	Double_t VxPull[3];
	TVectorD XvFit = Vfit->GetVtx();
	TMatrixDSym XvCov = Vfit->GetVtxCov();
	Int_t Ntr  = Vfit->GetNtrk();
	Double_t Chi2 = Vfit->GetVtxChi2();
	for(Int_t i=0; i<3; i++) VxPull[i] = (XvFit(i)-vTrue[i])/TMath::Sqrt(XvCov(i,i));
	h_Xv->Fill(VxPull[0]);
	h_Yv->Fill(VxPull[1]);
	h_Zv->Fill(VxPull[2]);
	//std::cout<<"Pulls: Vx = "<<VxPull[0]<<", "<<VxPull[1]<<", "<<VxPull[2]<<std::endl;
	// Chi/Ndof
	Double_t Ndof = 2 * (Double_t)Ntr;
	h_Chi2->Fill(Chi2/Ndof);
	//
	// Updated track parameter pulls
	Int_t Ntst = (Int_t) pTrue.size();
	if(Ntr != Ntst) {
		std::cout<<"VtxPulls:: mismatch on input track lists. Abort."<<std::endl;
		exit(-1);
	}else{
		Bool_t Mm = kTRUE;
		VertexMore* Vmore = new VertexMore(Vfit, Mm);
		// Track loop
		for(Int_t i=0; i<Ntr; i++){
			TVector3 pTrk(pTrue[i](0),pTrue[i](1),pTrue[i](2));
			TVector3 Xvm(vTrue[0]*1.e-3,vTrue[1]*1.e-3,vTrue[2]*1.e-3);
			Double_t Q = Vmore->GetCharge(i);
			TVectorD gPar = TrkUtil:: XPtoPar(Xvm, pTrk, Q, Bz);
			TVectorD GenPar = TrkUtil::ParToMm(gPar); // Back to mm
			TVectorD NewPar = Vfit->GetNewPar(i);
			TMatrixDSym ParCov = Vfit->GetNewCov(i);
			//
			// Track Parameter pulls
			Double_t ParPull[5];
			for(Int_t k=0; k<5; k++)ParPull[k] = (NewPar(k)-GenPar(k))/TMath::Sqrt(ParCov(k,k));
			h_D ->Fill(ParPull[0]);
			h_P0->Fill(ParPull[1]);
			h_C ->Fill(ParPull[2]);
			h_Z0->Fill(ParPull[3]);
			h_Ct->Fill(ParPull[4]);
			//
			// Get fit momenta momenta
			//
			TVector3 pRec = Vmore->GetMomentum(i);
			TMatrixDSym pCov = Vmore->GetMomentumC(i);
			Double_t MomPull[3];
			for(Int_t k=0; k<3; k++)MomPull[k] = (pRec(k)-pTrk(k))/TMath::Sqrt(pCov(k,k));
			h_Px->Fill(MomPull[0]);
			h_Py->Fill(MomPull[1]);
			h_Pz->Fill(MomPull[2]);
		}
		delete Vmore;
	}

}

// Print 
void VtxPulls::Print() 
{
	//
	// Vertex position pulls
	TCanvas * cnv = new TCanvas("cnv","Vertex position pulls",10,10, 600,600);
	cnv->Divide(2,2);
	// X
	cnv->cd(1); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_Xv->Fit("gaus"); h_Xv->Draw();
	// Y
	cnv->cd(2); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_Yv->Fit("gaus"); h_Yv->Draw();
	// Z
	cnv->cd(3); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_Zv->Fit("gaus"); h_Zv->Draw();
	// Chi2/Ndof
	cnv->cd(4); 
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_Chi2->Draw();
	//
	// Updated track parameter pulls
	TCanvas * cnv1 = new TCanvas("cnv1","Updated parameter pulls",50,50, 900,600);
	cnv1->Divide(3,2);
	// D
	cnv1->cd(1); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_D->Fit("gaus"); h_D->Draw();
	// Phi0
	cnv1->cd(2); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_P0->Fit("gaus"); h_P0->Draw();
	// C
	cnv1->cd(3); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_C->Fit("gaus"); h_C->Draw();
	// Z0
	cnv1->cd(4); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_Z0->Fit("gaus"); h_Z0->Draw();
	// Cotg Theta
	cnv1->cd(5); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_Ct->Fit("gaus"); h_Ct->Draw();
	//
	// Momenta pulls
	TCanvas * cnv2 = new TCanvas("cnv2","Momenta pulls",100,100, 900,600);
	cnv2->Divide(3,1);
	// Px
	cnv2->cd(1); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_Px->Fit("gaus"); h_Px->Draw();
	// Py
	cnv2->cd(2); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_Py->Fit("gaus"); h_Py->Draw();
	// Pz
	cnv2->cd(3); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_Pz->Fit("gaus"); h_Pz->Draw();
}
