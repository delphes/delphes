#include "B0KsKsPulls.h"
#include <TFile.h>
#include <TObjArray.h>
//
// Final state configuration description
//
// Constructor
B0KsKsPulls::B0KsKsPulls()
{
	// list
	Hlist(0);
	//
	// Pulls
	// 1st Ks vertex
	h_Ks1Xv = new TH1D("h_Ks1Xv","Ks1 X-vertex pull",100,-10.,10.);
	Hlist.Add(h_Ks1Xv);
	h_Ks1Yv = new TH1D("h_Ks1Yv","Ks1 Y-vertex pull",100,-10.,10.);
	Hlist.Add(h_Ks1Yv);
	h_Ks1Zv = new TH1D("h_Ks1Zv","Ks1 Z-vertex pull",100,-10.,10.);
	Hlist.Add(h_Ks1Zv);
	// 1st Ks vertex momentum
	h_Ks1Px = new TH1D("h_Ks1Px","Ks1 X-momentum pull",100,-10.,10.);
	Hlist.Add(h_Ks1Px);
	h_Ks1Py = new TH1D("h_Ks1Py","Ks1 Y-momentum pull",100,-10.,10.);
	Hlist.Add(h_Ks1Py);
	h_Ks1Pz = new TH1D("h_Ks1Pz","Ks1 Z-momentum pull",100,-10.,10.);
	Hlist.Add(h_Ks1Pz);
	// 2nd Ks vertex
	h_Ks2Xv = new TH1D("h_Ks2Xv","Ks2 X-vertex pull",100,-10.,10.);
	Hlist.Add(h_Ks2Xv);
	h_Ks2Yv = new TH1D("h_Ks2Yv","Ks2 Y-vertex pull",100,-10.,10.);
	Hlist.Add(h_Ks2Yv);
	h_Ks2Zv = new TH1D("h_Ks2Zv","Ks2 Z-vertex pull",100,-10.,10.);
	Hlist.Add(h_Ks2Zv);
	// 2nd Ks vertex momentum
	h_Ks2Px = new TH1D("h_Ks2Px","Ks2 X-momentum pull",100,-10.,10.);
	Hlist.Add(h_Ks2Px);
	h_Ks2Py = new TH1D("h_Ks2Py","Ks2 Y-momentum pull",100,-10.,10.);
	Hlist.Add(h_Ks2Py);
	h_Ks2Pz = new TH1D("h_Ks2Pz","Ks2 Z-momentum pull",100,-10.,10.);
	Hlist.Add(h_Ks2Pz);
	// B0 vertex
	h_B0Xv = new TH1D("h_B0Xv","B0 X-vertex pull",100,-10.,10.);
	Hlist.Add(h_B0Xv);
	h_B0Yv = new TH1D("h_B0Yv","B0 Y-vertex pull",100,-10.,10.);
	Hlist.Add(h_B0Yv);
	h_B0Zv = new TH1D("h_B0Zv","B0 Z-vertex pull",100,-10.,10.);
	Hlist.Add(h_B0Zv);
	// B0 vertex momentum
	h_B0Px = new TH1D("h_B0Px","B0 X-momentum pull",100,-10.,10.);
	Hlist.Add(h_B0Px);
	h_B0Py = new TH1D("h_B0Py","B0 Y-momentum pull",100,-10.,10.);
	Hlist.Add(h_B0Py);
	h_B0Pz = new TH1D("h_B0Pz","B0 Z-momentum pull",100,-10.,10.);
	Hlist.Add(h_B0Pz);

	//
	// Mass plots
	h_B0Mass  = new TH1D("h_B0Mass" , "B0 mass (GeV)" , 100, 5.1, 5.4);
	Hlist.Add(h_B0Mass);
	h_B0Mpull = new TH1D("h_B0Mpull", "B0 mass pull " , 100, -10., 10.);
	Hlist.Add(h_B0Mpull);
	h_B0Merr  = new TH1D("h_B0Merr" , "B0 mass error ", 100, 0., 0.15);
	Hlist.Add(h_B0Merr);

	// B0 flight path
	h_B0Lxy = new TH1D("h_B0Lxy", "B0 L 2D (fit - truth)", 100, -1.0, 1.0) ;
	h_B0LxyS = new TH1D("h_B0LxyS", "B0 L 2D error", 100, 0.0, 1.0) ;
	h_B0Lxyz = new TH1D("h_B0Lxyz", "B0 L 3D (fit - truth)", 100, -1.0, 1.0) ;
	h_B0LxyzS = new TH1D("h_B0LxyzS", "B0 L 3D error", 100, 0.0, 1.0) ;
	h_B0LxyPull = new TH1D("h_B0LxyPull", "B0 L 2D pull", 100, -10.0, 10.0) ;
	h_B0LxyzPull = new TH1D("h_B0LxyzPull", "B0 L 3D pull", 100, -10.0, 10.0) ;
}

//
// Destructor
B0KsKsPulls::~B0KsKsPulls() 
{
}

void B0KsKsPulls::Fill(VState* B0State, VertexMore* vKs1, VertexMore* vKs2, VertexMore* vB0)
{
	const Int_t nKs = 2;
	GenParticle* pB0    = B0State ->GetMother();	// Pointer to B0 genparticle
	VState* Ks1State    = B0State ->GetVState(0);	// 1st Ks
	GenParticle* pKs1   = Ks1State->GetMother();	// Pointer to Ks genparticle
	VState* Ks2State    = B0State ->GetVState(1);	// 2nd Ks
	GenParticle* pKs2   = Ks2State->GetMother();	// Pointer to Ks genparticle
	//
	// B0 Vertex
	TVector3 gvB0(pKs1->X, pKs1->Y, pKs1->Z);	// B0 generated vertex
	TVectorD rvB0    = vB0->GetXv();		// B0 reconstructed vertex
	TMatrixDSym cvB0 = vB0->GetXvCov();		// B0 vertex covariance
	//
	// B0 origin
	//
	TVector3 gvsB0(pB0->X, pB0->Y, pB0->Z);		// B0 start vertex
	TVector3 gLd = gvB0-gvsB0;			// B0 travel vector
	Double_t gLxy = gLd.Pt();	// Lxy generated
	Double_t gLxyz = gLd.Mag();	// Lxyz generated
	TVector3 rV(rvB0(0),rvB0(1), rvB0(2));
	TVector3 rLd = rV - gvsB0;	// Reconstructed travel vector
	Double_t rLxy = rLd.Pt();	// Reconstructed Lxy
	Double_t rLxyz = rLd.Mag();	// Reconstructed Lxyz
	// Flight path errors
	Double_t adLxy [2] = {gLd(0)/gLxy,  gLd(1)/gLxy};
	Double_t adLxyz[3] = {gLd(0)/gLxyz, gLd(1)/gLxyz, gLd(2)/gLxyz};
	TVectorD dLxy(2, adLxy);
	TVectorD dLxyz(3, adLxyz);
	//
	TMatrixDSym Bcov3 = cvB0;
	TMatrixDSym Bcov2 = cvB0.GetSub(0,1,0,1);
	h_B0Lxy->Fill(rLxy-gLxy);
	Double_t sLxy = TMath::Sqrt(Bcov2.Similarity(dLxy));
	h_B0LxyS->Fill(sLxy);
	h_B0LxyPull->Fill((rLxy-gLxy)/sLxy);
	h_B0Lxyz->Fill(rLxyz-gLxyz);
	Double_t sLxyz =  TMath::Sqrt(Bcov3.Similarity(dLxyz));
	h_B0LxyzS->Fill(sLxyz);
	h_B0LxyzPull->Fill((rLxyz-gLxyz)/sLxyz);
	// Fill
	h_B0Xv->Fill((rvB0(0)-gvB0(0))/TMath::Sqrt(cvB0(0,0)));
	h_B0Yv->Fill((rvB0(1)-gvB0(1))/TMath::Sqrt(cvB0(1,1)));
	h_B0Zv->Fill((rvB0(2)-gvB0(2))/TMath::Sqrt(cvB0(2,2)));
	// Bs Momentum
	TVector3 gpB0(pB0->Px, pB0->Py, pB0->Pz);	// B0 generated momentum
	TVector3 rpB0    = vB0->GetTotalP();		// B0 reconstructed momentum
	TMatrixDSym cpB0 = vB0->GetTotalPcov();		// B0 momentum covariance
	// Fill
	h_B0Px->Fill((rpB0(0)-gpB0(0))/TMath::Sqrt(cpB0(0,0)));
	h_B0Py->Fill((rpB0(1)-gpB0(1))/TMath::Sqrt(cpB0(1,1)));
	h_B0Pz->Fill((rpB0(2)-gpB0(2))/TMath::Sqrt(cpB0(2,2)));
	//
	// 1st Ks Vertex
	GenParticle* Pion1  = Ks1State->GetGen(0);	// Pi from first Ks
	TVector3 gvKs1(Pion1->X, Pion1->Y, Pion1->Z);	// K0s generated vertex
	TVectorD rvKs1 = vKs1->GetXv();			// Ks reconstructed vertex
	TMatrixDSym cvKs1 = vKs1->GetXvCov();		// Ks vertex covariance
	// Fill
	h_Ks1Xv->Fill((rvKs1(0)-gvKs1(0))/TMath::Sqrt(cvKs1(0,0)));
	h_Ks1Yv->Fill((rvKs1(1)-gvKs1(1))/TMath::Sqrt(cvKs1(1,1)));
	h_Ks1Zv->Fill((rvKs1(2)-gvKs1(2))/TMath::Sqrt(cvKs1(2,2)));
	// 1st Ks Momentum
	TVector3 gpKs1(pKs1->Px, pKs1->Py, pKs1->Pz);	// Ks generated momentum
	TVector3 rpKs1 = vKs1->GetTotalP();		// Ks reconstructed momentum
	TMatrixDSym cpKs1 = vKs1->GetTotalPcov();		// Ks momentum covariance
	// Fill
	h_Ks1Px->Fill((rpKs1(0)-gpKs1(0))/TMath::Sqrt(cpKs1(0,0)));
	h_Ks1Py->Fill((rpKs1(1)-gpKs1(1))/TMath::Sqrt(cpKs1(1,1)));
	h_Ks1Pz->Fill((rpKs1(2)-gpKs1(2))/TMath::Sqrt(cpKs1(2,2)));
	//
	// 2nd Ks Vertex
	GenParticle* Pion2  = Ks2State->GetGen(0);	// Pi from second Ks
	TVector3 gvKs2(Pion2->X, Pion2->Y, Pion2->Z);	// K0s generated vertex
	TVectorD rvKs2 = vKs2->GetXv();			// Ks reconstructed vertex
	TMatrixDSym cvKs2 = vKs2->GetXvCov();		// Ks vertex covariance
	// Fill
	h_Ks2Xv->Fill((rvKs2(0)-gvKs2(0))/TMath::Sqrt(cvKs2(0,0)));
	h_Ks2Yv->Fill((rvKs2(1)-gvKs2(1))/TMath::Sqrt(cvKs2(1,1)));
	h_Ks2Zv->Fill((rvKs2(2)-gvKs2(2))/TMath::Sqrt(cvKs2(2,2)));
	// 2nd Ks Momentum
	TVector3 gpKs2(pKs2->Px, pKs2->Py, pKs2->Pz);	// Ks generated momentum
	TVector3 rpKs2 = vKs2->GetTotalP();		// Ks reconstructed momentum
	TMatrixDSym cpKs2 = vKs2->GetTotalPcov();	// Ks momentum covariance
	// Fill
	h_Ks2Px->Fill((rpKs2(0)-gpKs2(0))/TMath::Sqrt(cpKs2(0,0)));
	h_Ks2Py->Fill((rpKs2(1)-gpKs2(1))/TMath::Sqrt(cpKs2(1,1)));
	h_Ks2Pz->Fill((rpKs2(2)-gpKs2(2))/TMath::Sqrt(cpKs2(2,2)));
	//
	// Bs mass plots
	TVector3 rB0Ks1 = vB0->GetMomentum(0);	// Momentum 1st Ks from B0
	TVector3 rB0Ks2 = vB0->GetMomentum(1);	// Momentum 2nd Ks from B0
	Double_t KsMass = pKs1->Mass;		// Ks mass
	Double_t B0Mass = pB0->Mass;		// B0 mass
	Double_t Ek1 = TMath::Sqrt(KsMass*KsMass+rB0Ks1.Mag2());	// Energy 1st Ks from B0
	Double_t Ek2 = TMath::Sqrt(KsMass*KsMass+rB0Ks2.Mag2());	// Energy 2nd Ks from B0
	Double_t Etot = Ek1 + Ek2;		// B0 energy
	TVector3 Ptot = rB0Ks1 + rB0Ks2;	// B0 momentum
	Double_t rB0Mass = TMath::Sqrt(Etot*Etot-Ptot.Mag2());	// Reconstructed B0 mass
	// Mass error
	TVector3 dm1 = (1./B0Mass)*((Etot/Ek1)*rB0Ks1-Ptot);	// dm/dp1
	TVector3 dm2 = (1./B0Mass)*((Etot/Ek2)*rB0Ks2-Ptot);	// dm/dp2
	Double_t dmp[6] = {dm1(0), dm1(1), dm1(2),dm2(0), dm2(1), dm2(2)};
	TVectorD DmDp(6, dmp);
	TMatrixDSym Bcov = vB0->GetBigCov();			
	TMatrixDSym pcov = Bcov.GetSub(3,8,3,8);		// cov(p1,p2)
	Double_t mB0Err = TMath::Sqrt(pcov.Similarity(DmDp));	// Error on B0 mass
	// Fill
	h_B0Mass->Fill(rB0Mass);
	h_B0Merr->Fill(mB0Err);
	h_B0Mpull->Fill((rB0Mass-B0Mass)/mB0Err);
}

// Print 
void B0KsKsPulls::Print() 
{

	//
	// 1st Ks plots
	TCanvas * cnv = new TCanvas("cnv","1st Ks vertex pulls",10,10, 900,600);
	cnv->Divide(3,2);
	// X
	cnv->cd(1); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_Ks1Xv->Fit("gaus"); h_Ks1Xv->Draw();
	// Y
	cnv->cd(2); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_Ks1Yv->Fit("gaus"); h_Ks1Yv->Draw();
	// Z
	cnv->cd(3); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_Ks1Zv->Fit("gaus"); h_Ks1Zv->Draw();
	// Px
	cnv->cd(4); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_Ks1Px->Fit("gaus"); h_Ks1Px->Draw();
	// Py
	cnv->cd(5); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_Ks1Py->Fit("gaus"); h_Ks1Py->Draw();
	// Pz
	cnv->cd(6); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_Ks1Pz->Fit("gaus"); h_Ks1Pz->Draw();
	//
	// 2nd  Ks plots
	TCanvas * cnv2 = new TCanvas("cnv2","2nd Ks vertex pulls",50,50, 900,600);
	cnv2->Divide(3,2);
	// X
	cnv2->cd(1); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_Ks2Xv->Fit("gaus"); h_Ks2Xv->Draw();
	// Y
	cnv2->cd(2); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_Ks2Yv->Fit("gaus"); h_Ks2Yv->Draw();
	// Z
	cnv2->cd(3); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_Ks2Zv->Fit("gaus"); h_Ks2Zv->Draw();
	// Px
	cnv2->cd(4); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_Ks2Px->Fit("gaus"); h_Ks2Px->Draw();
	// Py
	cnv2->cd(5); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_Ks2Py->Fit("gaus"); h_Ks2Py->Draw();
	// Pz
	cnv2->cd(6); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_Ks2Pz->Fit("gaus"); h_Ks2Pz->Draw();
	//
	// B0 plots
	TCanvas * cnv1 = new TCanvas("cnv1","B0 vertex pulls",100,100, 900,600);
	cnv1->Divide(3,2);
	// X
	cnv1->cd(1); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_B0Xv->Fit("gaus"); h_B0Xv->Draw();
	// Y
	cnv1->cd(2); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_B0Yv->Fit("gaus"); h_B0Yv->Draw();
	// Z
	cnv1->cd(3); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_B0Zv->Fit("gaus"); h_B0Zv->Draw();
	// Px
	cnv1->cd(4); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_B0Px->Fit("gaus"); h_B0Px->Draw();
	// Py
	cnv1->cd(5); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_B0Py->Fit("gaus"); h_B0Py->Draw();
	// Pz
	cnv1->cd(6); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_B0Pz->Fit("gaus"); h_B0Pz->Draw();
	//
	// B0 mass plots
	TCanvas * cnv3 = new TCanvas("cnv3","B0 mass plots",200,200, 900,600);
	cnv3->Divide(2,2);
	cnv3->cd(1); 
	gStyle->SetOptStat(111111); 
	h_B0Mass->Draw();
	cnv3->cd(2); 
	gStyle->SetOptStat(111111); 
	h_B0Merr->Draw();
	cnv3->cd(3); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_B0Mpull->Fit("gaus"); h_B0Mpull->Draw();
	//
	// Flight path plots
	TCanvas * cnvp = new TCanvas("cnvp","B0 flight path plots",250,250, 900,600);
	cnvp->Divide(3,2);
	cnvp->cd(1);
	h_B0Lxy->Draw();
	cnvp->cd(2);
	h_B0LxyS->Draw(); gPad->SetLogy(1);
	cnvp->cd(3); gPad->SetLogy(1); gStyle->SetOptFit(1111);
	h_B0LxyPull->Fit("gaus"); h_B0LxyPull->Draw(); 
	cnvp->cd(4);
	h_B0Lxyz->Draw();
	cnvp->cd(5);
	h_B0LxyzS->Draw(); gPad->SetLogy(1);
	cnvp->cd(6); gPad->SetLogy(1);gStyle->SetOptFit(1111);
	h_B0LxyzPull->Fit("gaus"); h_B0LxyzPull->Draw(); 
	//
	// write out everything
	TFile hFile("hB0KsKs.root","RECREATE");
	Hlist.Write();
	hFile.Close();
}
