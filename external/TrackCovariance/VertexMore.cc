#include "TrkUtil.h"
#include <iostream>
#include "VertexMore.h"

// Constructors
VertexMore::VertexMore(VertexFit* V)
{
	fUnits = kFALSE;
	fV = V;
	fNtr = fV->GetNtrk();
	fPar.ResizeTo(5);	fPar.Zero();
	fCov.ResizeTo(5, 5);	fCov.Zero();
	fCp.ResizeTo(3, 3);
	CalcParCov();
	fXv.ResizeTo(3);
	fXv = fV->GetVtx();
	if(fQtot != 0){
		fPar = MakeVpar();
		fCov = MakeVcov();
	}
	fBigCov.ResizeTo(3 * (fNtr + 1), 3 * (fNtr + 1));
	FillBigCov();
}
VertexMore::VertexMore(VertexFit* V, Bool_t opt)
{
	fUnits = opt;
	fV = V;
	fNtr = fV->GetNtrk();
	fPar.ResizeTo(5);	fPar.Zero();
	fCov.ResizeTo(5, 5);	fCov.Zero();
	fCp.ResizeTo(3, 3);
	CalcParCov();
	fXv.ResizeTo(3);
	fXv = fV->GetVtx();
	if(fQtot != 0){
		fPar = MakeVpar();
		fCov = MakeVcov();
	}
	fBigCov.ResizeTo(3 * (fNtr + 1), 3 * (fNtr + 1));
	FillBigCov();
}
//
// Destructor
VertexMore::~VertexMore()
{
	for (Int_t i = 0; i < fNtr; i++) {
		fpi[i]->Clear();	
		fCpi[i]->Clear();	
	}
	fpi.clear();
	fCpi.clear();
	fQ.clear();
}
//
// Main method to calculate parameters and covariance from vertex

TMatrixD VertexMore::dPdAlf(Int_t i)
{
	TMatrixD dPdPar(5, 3); dPdPar.Zero();
	TVector3 p = GetMomentum(i); 
	TVectorD par = fV->GetNewPar(i);
	Double_t C = par(2);
	//
	dPdPar(1, 0) = -p(1);
	dPdPar(1, 1) = p(0);
	dPdPar(2, 0) = -p(0) / C;
	dPdPar(2, 1) = -p(1) / C;
	dPdPar(2, 2) = -p(2) / C;
	dPdPar(4, 2) = p.Pt();
	//
	return dPdPar;
}
void VertexMore::CalcParCov()
{
	//
	// Set magnetic field and units
	Double_t Bz = 2.;
	SetB(Bz);
	fa = -Bz * cSpeed();
	if(fUnits)fa = -Bz * cSpeed()*1.0e-3;
	//
	// Vertex position
	TVectorD xx = fV->GetVtx();
	TVector3 xv(xx.GetMatrixArray());
	//
	// Vertex momentum
	TVector3 pv(0., 0., 0.);
	//
	// Vertex charge
	fQtot = 0.;
	//
	// Vertex momentum error
	TMatrixD CovPtot(3, 3); CovPtot.Zero();
	std::vector<TMatrixD*> dPdAlfa;
	//
	for (Int_t i = 0; i < fNtr; i++) {
		// Unpack ith track parameters
		TVectorD par = fV->GetNewPar(i);
		//
		Double_t C = par(2);
		Double_t Q = TMath::Sign(1., -C * Bz);
		fQ.push_back(Q);
		fQtot += Q;
		//
		Double_t D = par(0);
		Double_t ph0 = par(1);
		Double_t z_0 = par(3);
		Double_t lm = par(4);
		Double_t a = fa * fQ[i];
		//
		// Store track momenta
		Double_t pt = a / (2. * C);
		Double_t s = 2 * TMath::ASin(C * TMath::Sqrt(TMath::Max(xv.Perp2() - D * D, 0.) / (1. + 2. * C * D)));
		TVector3 p(pt * TMath::Cos(s + ph0), pt * TMath::Sin(s + ph0), pt*lm);
		//
		fpi.push_back(new TVector3(p));
		//
		pv += p;	// Total momentum
		//
		// Get momentum errors
		TMatrixD dPdPar = dPdAlf(i);
		dPdAlfa.push_back(new TMatrixD(dPdPar));
		TMatrixDSym ParCov = fV->GetNewCov(i);
		TMatrixDSym pCov = ParCov.SimilarityT(dPdPar);
		fCpi.push_back(new TMatrixDSym(pCov));
		//
	}
	//
	// Total momentum
	fP = pv;
	//
	// Total momentum error
	for (Int_t i = 0; i < fNtr; i++) {
		for (Int_t j = 0; j < fNtr; j++) {
			TMatrixD Cij = fV->GetNewCov(i, j);
			TMatrixD dPi = *dPdAlfa[i];
			TMatrixD dPj = *dPdAlfa[j];
			TMatrixD dPiT(TMatrixD::kTransposed, dPi);
			CovPtot += dPiT * (Cij * dPj);
		}
	}
	//
	// Symmetrize Total momentum matrix
	for (Int_t i = 0; i < 3; i++) {
		for (Int_t j = 0; j < 3; j++)fCp(i, j) = 0.5 * (CovPtot(i, j) + CovPtot(j, i));
	}
	// Clean
	for (Int_t i = 0; i < fNtr; i++) {
		dPdAlfa[i]->Clear();
		delete dPdAlfa[i];
	}
	dPdAlfa.clear();
	//
}
//
// Vertex parameters
TVectorD VertexMore::MakeVpar()
{
	//
	// Vertex position
	TVectorD xx = fV->GetVtx();
	TVector3 xv(xx.GetMatrixArray());
	if(fUnits) xv *= 1.0e-3;	// Change to meters
	TVectorD Par(5); Par.Zero();
	if(fQtot != 0.0)Par = XPtoPar(xv, fP, fQtot);
	else {
		std::cout << "VertexMore::GetVpar: zero charge vertex not supported" << std::endl;
		std::exit(1);
	}
	fPar = Par;
	if(fUnits)fPar = ParToMm(Par);	// Back to mm
	//
	return fPar;
}
//
// Vertex parameter errors
//
TMatrixD VertexMore::DparDx(TVector3 xv, TVector3 pv, Double_t Q)
{
	//
	// Derivative of track parameters wrt point on  track
	//
	if (Q == 0.0) {
		std::cout << "VertexMore::DparDx: zero charge vertex not supported" << std::endl;
		std::exit(1);
	}
	// Track parameters
	TVectorD Par = XPtoPar(xv, pv, Q);
	Double_t D = Par(0);
	Double_t ph0 = Par(1);
	Double_t C = Par(2);
	Double_t z_0 = Par(3);
	Double_t lm = Par(4);
	Double_t a = fa * Q;
	//
	// Derivative matrix
	TMatrixD dParX(5,3); dParX.Zero();
	//
	// dT/x, dT/p
	Double_t pt = pv.Pt();
	Double_t R2 = xv.Pt() * xv.Pt();
	Double_t T = TMath::Sqrt(pt * pt - 2. * a * (xv.X() * pv.Y() - xv.Y() * pv.X()) + a * a * R2);
	TVectorD dTdx(3);
	dTdx(0) = a * (-pv.Y() + a * xv.X()) / T;
	dTdx(1) = a * (pv.X() + a * xv.Y()) / T;
	dTdx(2) = 0.;
	// D derivatives
	for (Int_t i = 0; i < 2; i++)dParX(0, i) = dTdx(i) / a;
	// Phi0 derivatives
	Double_t tgp = TMath::Tan(ph0);
	Double_t cs2 = pow(TMath::Cos(ph0), 2);
	dParX(1, 0) = -(a / (pv.X() + a * xv.Y())) * cs2;
	dParX(1, 1) = -(a / (pv.X() + a * xv.Y())) * cs2 * tgp;
	// Z_0 derivatives
	// Use approximation z0 ~ z - lm*sqrt(R^2-D^2)
	dParX(3, 0) = -(lm / TMath::Sqrt(R2 - D * D)) * (xv.X() - D * dParX(0, 0));
	dParX(3, 1) = -(lm / TMath::Sqrt(R2 - D * D)) * (xv.Y() - D * dParX(0, 1));
	dParX(3, 2) = 1.0;
	//
	return dParX;
}
//
void VertexMore::FillBigCov()
{
	//
	// Fill x vertex track momenta correlations
	//
	TMatrixDSub(fBigCov, 0, 2, 0, 2) = fV->GetVtxCov(); // <x.x>
	TMatrixD CovPX(3, 3); CovPX.Zero();
	for (Int_t i = 0; i < fNtr; i++) {
		TMatrixD dpda = dPdAlf(i);	// (5,3)
		for (Int_t k = 0; k < fNtr; k++) {
			TMatrixD dada0 = fV->DaiDa0k(i, k); // (5,5)
			TMatrixDSym Cv = fV->GetOldCov(k);  // (5,5)
			TMatrixD dxda0 = fV->GetDxvDpar0(k);// (3,5)
			TMatrixD dxda0t(TMatrixD::kTransposed, dxda0); // (5,3)
			TMatrixD dpdat(TMatrixD::kTransposed, dpda);
			CovPX += (dpdat * dada0) * (Cv * dxda0t); // (3.5)(5,5)(5,5)(5,3)
		}
		TMatrixD CovXP(TMatrixD::kTransposed, CovPX);
		TMatrixDSub(fBigCov, 3 * (i + 1), 3 * (i + 2) - 1, 0, 2) = CovXP; // <x.p>
		TMatrixDSub(fBigCov, 0, 2, 3 * (i + 1), 3 * (i + 2) - 1) = CovPX; // <p.x>
	}
	//
	// Fill momenta correlations
	for (Int_t i = 0; i < fNtr; i++) {
		for (Int_t j = 0; j < fNtr; j++) {
			if(i == j)TMatrixDSub(fBigCov, 3 * (i + 1), 3 * (i + 2) - 1, 3 * (i + 1), 3 * (i + 2) - 1) = GetMomentumC(i); // <pi.pi>
			else {
				TMatrixD dPdai = dPdAlf(i); // (5,3)
				TMatrixD dPdait(TMatrixD::kTransposed, dPdai); // (3,5)
				TMatrixD dPdaj = dPdAlf(j); // (5,3)
				TMatrixD Cov_ij = fV->GetNewCov(i, j); // (5,5)
				TMatrixD Pij = dPdait * (Cov_ij * dPdaj); // (3,3)
				TMatrixDSub(fBigCov, 3 * (i + 1), 3 * (i + 2) - 1, 3 * (j + 1), 3 * (j + 2) - 1) = Pij;
			}
		}
	}

}
//
TMatrixD VertexMore::DparDp(TVector3 xv, TVector3 pv, Double_t Q)
{
	//
	// Derivative of track parameters wrt point on track
	//
	if (Q == 0.0) {
		std::cout << "VertexMore::DparDp: zero charge vertex not supported" << std::endl;
		std::exit(1);
	}
	// Track parameters
	TVectorD Par = XPtoPar(xv, pv, Q);
	Double_t D = Par(0);
	Double_t ph0 = Par(1);
	Double_t C = Par(2);
	Double_t z_0 = Par(3);
	Double_t lm = Par(4);
	Double_t a = fa * Q;
	//
	// Derivative matrix
	TMatrixD dParP(5, 3); dParP.Zero();
	//
	// dT/x, dT/p
	Double_t pt = pv.Pt();
	Double_t R2 = xv.Pt() * xv.Pt();
	Double_t T = TMath::Sqrt(pt * pt - 2. * a * (xv.X() * pv.Y() - xv.Y() * pv.X()) + a * a * R2);
	TVectorD dTdp(3);
	dTdp(0) = (pv.X() + a * xv.Y()) / T;
	dTdp(1) = (pv.Y() - a * xv.X()) / T;
	dTdp(2) = 0.;
	// D derivatives
	for (Int_t i = 0; i < 2; i++)dParP(0, i) = (dTdp(i) - pv(i) / pt) / a;
	// Phi0 derivatives
	Double_t tgp = TMath::Tan(ph0);
	Double_t cs2 = pow(TMath::Cos(ph0), 2);
	dParP(1, 0) = -tgp * cs2 / (pv.X() + a * xv.Y());
	dParP(1, 1) = cs2 / (pv.X() + a * xv.Y());
	// C derivatives
	dParP(2, 0) = -a * pv.X() / (2 * pt * pt * pt);
	dParP(2, 1) = -a * pv.Y() / (2 * pt * pt * pt);
	// lambda derivatives
	dParP(4, 0) = -pv.Z() * pv.X() / (pt * pt * pt);
	dParP(4, 1) = -pv.Z() * pv.Y() / (pt * pt * pt);
	dParP(4, 2) = 1.0 / pt;
	// Z_0 derivatives
	// Use approximation z0 ~ z - lm*sqrt(R^2-D^2)
	dParP(3, 0) = -dParP(4, 0) * TMath::Sqrt(R2 - D * D) + (lm / TMath::Sqrt(R2 - D * D)) * D * dParP(0, 0);
	dParP(3, 1) = -dParP(4, 1) * TMath::Sqrt(R2 - D * D) + (lm / TMath::Sqrt(R2 - D * D)) * D * dParP(0, 1);
	dParP(3, 2) = -dParP(4, 2) * TMath::Sqrt(R2 - D * D);
	//
	return dParP;
}
//

TMatrixDSym VertexMore::MakeVcov()
{
	// da = da/dx dx + da/dp dp
	// <da.da'> = da/dx<x.x'>(da/dx)'+da/dx<x.p'>(da/dp)'+
	//	      da/dp<p.p'>(da/dp)'+da/dp<p.x'>(da/dx)'
	//
	if (fQtot == 0.0){
		std::cout << "VertexMore::GetVcov: zero charge vertex not supported" << std::endl;
		std::exit(1);
	}
	//Double_t a = fa*fQtot;
	//
	// Vertex position
	TVectorD xx = fV->GetVtx();
	TVector3 xv(xx.GetMatrixArray());
	TVector3 pv = GetTotalP();
	//
	// Vertex parameters
	Double_t D   = fPar(0);
	Double_t ph0 = fPar(1);
	Double_t C   = fPar(2);
	Double_t z_0 = fPar(3);
	Double_t lm  = fPar(4);
	//
	// Parameter derivatives dPar/dx, dPar/dp
	TMatrixD dParX = DparDx(xv, pv, fQtot);
	TMatrixD dParP = DparDp(xv, pv, fQtot);
	//
	Double_t s = 2 * TMath::ASin(C * TMath::Sqrt(TMath::Max(xv.Perp2() - D * D, 0.) / (1. + 2. * C * D)));
	TMatrixD dXdPar = derXdPar(fPar, s);
	// 
	// <p'.x> = {S_ik dp_i/da_i da_i/da0_k C_k At_iD_i}D^-1
	TMatrixD CovPX(3, 3); CovPX.Zero();
	//
	for (Int_t i = 0; i < fNtr; i++) {
		TMatrixD dpda = dPdAlf(i);	// (5,3)
		for (Int_t k = 0; k < fNtr; k++) {
			TMatrixD dada0 = fV->DaiDa0k(i, k); // (5,5)
			TMatrixDSym Cv = fV->GetOldCov(k);  // (5,5)
			TMatrixD dxda0 = fV->GetDxvDpar0(k);// (3,5)
			TMatrixD dxda0t(TMatrixD::kTransposed, dxda0); // (5,3)
			TMatrixD dpdat(TMatrixD::kTransposed,dpda);
			CovPX += (dpdat * dada0) * (Cv * dxda0t); // (3.5)(5,5)(5,5)(5,3)
		}
	}
	//
	TMatrixDSym Xcov = fV->GetVtxCov();
	fCov += Xcov.Similarity(dParX);
	TMatrixDSym Pcov = GetTotalPcov();
	fCov += Pcov.Similarity(dParP);
	TMatrixD dParXt(TMatrixD::kTransposed, dParX);
	TMatrixD XPcross = dParP * (CovPX * dParXt);
	TMatrixDSym XPsym(5);
	for (Int_t i = 0; i < 5; i++) {
		for (Int_t j = 0; j < 5; j++)XPsym(i, j) = XPcross(i, j) + XPcross(j, i);
	}
	fCov += XPsym;
	//
	return fCov;
}
//

