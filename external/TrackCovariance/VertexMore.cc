#include "VertexMore.h"

void VertexMore::Mprt(TMatrixD M, Bool_t Opt)
{
	Int_t N = M.GetNrows();
	TMatrixD Mn = M;
	TVectorD d(N);
	if(Opt){			// Normalize if chosen
		for(Int_t i=0; i<N; i++){
			if(M(i,i)<= 0.0){
				std::cout<<"M("<<i<<", "<<i<<") = "<<M(i,i)<<std::endl;
				d(i) = TMath::Abs(M(i,i)) + 1.e-6;
			} else 	d(i) = TMath::Sqrt(M(i,i));
		}
		for(Int_t i=0; i<N; i++){
			for(Int_t j=0; j<N; j++)Mn(i,j) = M(i,j)/(d(i)*d(j));
		}
	}
	std::cout<<std::endl;
	for(Int_t i = 0; i<N; i++){
		for(Int_t j = 0; j<N; j++)std::cout<<Mn(i,j)<<"\t";
		std::cout<<std::endl;
	}
}
// Constructors
//
void VertexMore::Init(VertexFit* V)
{
	//
	// Set magnetic field and units
	Double_t Bz = 2.;
	SetB(Bz);
	fa = -Bz * cSpeed();
	if(fUnits)fa = -Bz * cSpeed()*1.0e-3;
	//
	// Initialize various arrays
	//
	fV = V;
	fNtr = fV->GetNtrk();
	fXv.ResizeTo(3);
	fXvCov.ResizeTo(3,3);
	fXv = fV->GetVtx();
	fXvCov = fV->GetVtxCov();
	fPar.ResizeTo(5);	fPar.Zero();
	fCov.ResizeTo(5, 5);	fCov.Zero();
	fCp.ResizeTo(3, 3);
	CalcParCov();
	fBigCov.ResizeTo(3 * (fNtr + 1), 3 * (fNtr + 1));
	FillBigCov();
	fBigPar.ResizeTo(3 * (fNtr + 1));
	FillBigPar();
	fPar = MakeVpar();
	fCov = MakeVcov();
	fNc = 0;
}
//
VertexMore::VertexMore(VertexFit* V)
{
	fUnits = kFALSE;
	Init(V);
}
VertexMore::VertexMore(VertexFit* V, Bool_t opt)
{
	fUnits = opt;
	Init(V);
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
//
Double_t VertexMore::dSdD(Int_t i)	// ***** NOT USED *****
{
	// Derivative of phase of ith track wrt D
	//
	TVector3 p = GetMomentum(i); 
	TVectorD xv = fV->GetVtx();
	Double_t R2 = xv(0)*xv(0)+xv(1)*xv(1);
	TVectorD par = fV->GetNewPar(i);
	Double_t D = par(0);
	Double_t C = par(2);
	//
	Double_t A = C*TMath::Sqrt(TMath::Max(0.0,R2-D*D)/(1.+2.*C*D));
	Double_t Num = -2.*C*C*(D+C*(R2+D*D));
	Double_t Den = (1.+2.*C*D)*(1.+2.*C*D)*A*TMath::Sqrt(1.-A*A);
	return Num/Den;
}

Double_t VertexMore::dSdC(Int_t i)	// ***** NOT USED *****
{
	// Derivative of phase of ith track wrt C
	//
	//TVector3 p = GetMomentum(i); 
	TVectorD xv = fV->GetVtx();
	Double_t R2 = xv(0)*xv(0)+xv(1)*xv(1);
	TVectorD par = fV->GetNewPar(i);
	Double_t D = par(0);
	Double_t C = par(2);
	//
	Double_t A = C*TMath::Sqrt(TMath::Max(0.0,R2-D*D)/(1.+2.*C*D));
	Double_t Num = 2.*A*(1.+C*D);
	Double_t Den = C*(1.+2.*C*D)*TMath::Sqrt(1.-A*A);
	return Num/Den;
}

TMatrixD VertexMore::dPdX(Int_t i)	// ***** NOT USED *****
{
	TVectorD par = fV->GetNewPar(i);
	TVectorD xv = fV->GetVtx();
	//
	Double_t C = par(2);
	Double_t z_0 = par(3);
	Double_t lm = par(4);
	//
	TVector3 p = GetMomentum(i);
	//
	TVectorD sx = dsdx(xv, par);
	//
	TMatrixD dPX(3,3); dPX.Zero();
	// px
	dPX(0,0) = -p(1)*sx(0);	// x
	dPX(1,0) = -p(1)*sx(1);	// y
	// py
	dPX(0,1) =  p(0)*sx(0);	// x
	dPX(1,1) =  p(0)*sx(1);	// y
	//
	return dPX;
}

TMatrixD VertexMore::dPdAlf(Int_t i)
{
	TMatrixD dPdPar(5, 3); dPdPar.Zero();
	TVectorD par = fV->GetNewPar(i);
	TVectorD xv = fV->GetVtx();
	Double_t ph0 = par(1);
	Double_t z_0 = par(3);
	Double_t lm = par(4);
	
	TVector3 p = GetMomentum(i); 

	if(fV->IsCharged(i)){		//Charged
		Double_t C = par(2);
		//
		// px
		dPdPar(1, 0) = -p(1);			// phi0
		dPdPar(2, 0) = -p(0)/C;			// C
		// py
		dPdPar(1, 1) =  p(0);			// Phi0
		dPdPar(2, 1) = -p(1)/C;			// C
		// pz
		dPdPar(2, 2) = -p(2) / C;		// C
		dPdPar(4, 2) = p.Pt();			// cot(theta)
	}else{				// Neutral
		Double_t pt = p.Pt();
		// 
		// px
		dPdPar(1, 0) = -p(1);		// phi0		
		dPdPar(2, 0) = TMath::Cos(ph0);	// pt
		// py
		dPdPar(1, 1) =  p(0);		// phi0
		dPdPar(2, 1) = TMath::Sin(ph0);	// pt
		// pz
		dPdPar(2, 2) = lm;		// pt
		dPdPar(4, 2) = pt;		// cot(theta)
	}
	//
	return dPdPar;
}

TVectorD VertexMore::dPds(Int_t i)
{
	TVectorD dPs(3); dPs.Zero();
	TVector3 p = GetMomentum(i); 
	//
	if(fV->IsCharged(i)){				//Charged
		dPs(0) = -p(1);
		dPs(1) =  p(0);
		dPs(2) =    0.;
	}
	//
	return dPs;
}
//
TMatrixD VertexMore::dXdAlf(Int_t i)		// *** NOT USED ***
{
	TMatrixD dXdPar(5, 3); dXdPar.Zero();
	TVectorD xv = fV->GetVtx();
	Double_t R2 = xv(0)*xv(0)+xv(1)*xv(1);
	TVectorD par = fV->GetNewPar(i);
	Double_t D = par(0);
	Double_t ph = par(1);
	Double_t sf = TMath::Sin(ph);
	Double_t cf = TMath::Cos(ph);
	Double_t z0 = par(3);
	Double_t ct = par(4);

	if(fV->IsCharged(i)){			// Charged
		//
		Double_t C = par(2);
		Double_t s = GetPhase(xv, par);
		//
		TVectorD dsdpr = dsdPar(xv, par);
		Double_t sd = dsdpr(0);
		Double_t sp0= dsdpr(1);
		Double_t sc = dsdpr(2);
		//
		Double_t sfs= TMath::Sin(s+ph);
		Double_t cfs= TMath::Cos(s+ph);
		//
		// x
		dXdPar(0, 0) = -sf+cfs*sd;		// D
		dXdPar(1, 0) = -D*cf+(cfs-cf)/(2.* C);	// phi0
		dXdPar(2, 0) = -(sfs-sf)/(2.*C*C)+cfs/(2.*C)*sc;	// C
		// y
		dXdPar(0, 1) =  cf+sfs*sd;		// D
		dXdPar(1, 1) = -D*sf+(sfs-sf)/(2.*C);	// phi0
		dXdPar(2, 1) =  (cfs-cf)/(2.*C*C)+sfs/(2.*C)*sc ;	// C
		// z
		dXdPar(2, 2) = -ct*s/(2.*C*C)+ct*s*sc/(2.*C);		// C
		dXdPar(3, 2) = 	1.;			// z0
		dXdPar(4, 2) =  s/(2.*C);		// ctg
	}
	else{					// Neutral
		Double_t s = xv(0)*cf+xv(1)*sf;
		//
		// x
		dXdPar(0, 0) = -sf;			// D
		//dXdPar(1, 0) = -D*cf - s*sf;		// phi0
		dXdPar(1, 0) = -s*sf;			// phi0
		// y
		dXdPar(0, 1) =  cf;			// D
		//dXdPar(1, 1) = -D*sf + s*cf;		// phi0
		dXdPar(1, 1) =  s*cf;			// phi0
		// z
		dXdPar(1, 2) =  D*ct;			// phi0
		dXdPar(3, 2) = 	1.;			// z0
		dXdPar(4, 2) =  s;			// ctg	
	}
	//
	return dXdPar;
}
//
TMatrixD VertexMore::DpDa0(Int_t i, Int_t k)
{
	//
	//	Momentum derivatives wrt initial track parameters
	//
	TMatrixD dPidPar0k(3,5);		// Return matrix
	// Get momentum derivatives wrt track parameters
	TMatrixD dPdPar_i = dPdAlf(i);		// (5, 3)
	TMatrixD dPdPar_iT(TMatrixD::kTransposed, dPdPar_i);	// (3, 5)
	// Final track parameter derivative wrt initial parameters
	TMatrixD DaDa0_k = fV->DaiDa0k(i, k);	// (5,5)
	// 
	// Update return matrix
	dPidPar0k = dPdPar_iT*DaDa0_k;		// Track par. components
	//
	// Get derivatives wrt phase
	TVectorD dPdph_i = dPds(i);		// (3)
	// Phase derivative wrt initial track parameters
	TVectorD DsDa0_k = fV->DsiDa0k(i,k);	// (5)
	//
	// Update return matrix
	TMatrixD sPart(3,5);
	dPidPar0k += sPart.Rank1Update(dPdph_i, DsDa0_k,1.0);
	//
	return dPidPar0k;
}
//
TMatrixD VertexMore::BuildPPcov(Int_t i, Int_t j)
{
	TMatrixD pCov(3,3); pCov.Zero();
	//
	// Loop over all initial track parameters
	for(Int_t k=0; k<fNtr; k++){
		TMatrixD Dpi  = DpDa0(i, k);
		TMatrixD Dpj  = DpDa0(j, k);
		TMatrixD DpjT(TMatrixD::kTransposed, Dpj);
		TMatrixDSym Ck = fV->GetOldCov(k);	// Initial track par. cov.
		pCov += Dpi*(Ck*DpjT);
	}
	return pCov;
}
//
TMatrixD VertexMore::BuildPXcov(Int_t i)
{
	TMatrixD pxCov(3,3); pxCov.Zero();
	//
	// Loop over all initial track parameters
	for(Int_t k=0; k<fNtr; k++){
		TMatrixD Dpi  = DpDa0(i, k);
		TMatrixDSym Ck = fV->GetOldCov(k);	// Initial track par. cov.
		TMatrixD DxvDa0_k = fV->GetDxvDpar0(k);
		TMatrixD DxvDa0_kT(TMatrixD::kTransposed, DxvDa0_k);
		pxCov += Dpi*(Ck*DxvDa0_kT);
	}
//
	return pxCov;
}
//
void VertexMore::CalcParCov()
{
	//
	TVectorD xv = fXv;		// Vertex position
	TVector3 pv(0., 0., 0.);	// Vertex momentum
	fQtot = 0.;			// Vertex charge
	//
	// Vertex momentum error
	//
	for (Int_t i = 0; i < fNtr; i++) {
		// Unpack ith track parameters
		TVectorD par = fV->GetNewPar(i);
		//Double_t D = par(0);
		Double_t ph0 = par(1);
		Double_t z_0 = par(3);
		Double_t lm = par(4);
		//
		Double_t Q = 0.;
		TVector3 p;
		if(fV->IsCharged(i)){			// Charged
			Double_t C = par(2);
			Q = TMath::Sign(1., C * fa);
			Double_t a = fa * Q;
			//
			// Store track momenta
			Double_t pt = TMath::Abs(a / (2. * C));
			Double_t s = fV->GetPhase(i);
			TVector3 pc(pt * TMath::Cos(s + ph0), pt * TMath::Sin(s + ph0), pt*lm);
			p = pc;
		}else{					// Neutral
			Double_t pt = par(2);
			TVector3 pc(pt * TMath::Cos(ph0), pt * TMath::Sin(ph0), pt*lm);
			p = pc;
		}
		fQ.push_back(Q);
		fQtot += Q;
		//
		fpi.push_back(new TVector3(p));
		//
		pv += p;	// Total momentum
		//
		// Get momentum errors
		//
		TMatrixDSym pCov(3);
		TMatrixD pCovTmp = BuildPPcov(i, i);	// Get momentum covariance
		for(Int_t h1=0; h1<3; h1++){		// Symmetrize
			for(Int_t h2=0; h2<3; h2++)pCov(h1,h2) = 0.5*(pCovTmp(h1,h2)+pCovTmp(h2,h1));
		}
		fCpi.push_back(new TMatrixDSym(pCov));
		//
	}
	//
	// Total momentum
	fP = pv;
	//
	// Total momentum error
	TMatrixD CovPtot(3, 3); CovPtot.Zero();
	for (Int_t i = 0; i < fNtr; i++) {
		for (Int_t j = 0; j < fNtr; j++) {
			CovPtot += BuildPPcov(i,j);
		}
	}
	//
	// Symmetrize Total momentum matrix
	for (Int_t i = 0; i < 3; i++) {
		for (Int_t j = 0; j < 3; j++)fCp(i, j) = 0.5 * (CovPtot(i, j) + CovPtot(j, i));
	}
	//
}
//
// Vertex parameters
TVectorD VertexMore::MakeVpar()
{
	//
	// Vertex position
	TVector3 xv(fXv(0), fXv(1), fXv(2));
	TVectorD Par(5); Par.Zero();
	if(fQtot != 0.0){			// Charged
		if(fUnits) xv *= 1.0e-3;	// Change to meters
		Par = XPtoPar(xv, fP, fQtot);
		if(fUnits)Par = ParToMm(Par);	// Back to mm
	}	
	else Par = XPtoPar_N(xv, fP);		// Neutral
	fPar = Par;
	//
	return fPar;
}
//
void VertexMore::FillBigCov()
{
	//
	// Fill x vertex track momenta correlations
	//
	Int_t Nbig = 3*(fNtr+1);
	TMatrixD BigCov(Nbig,Nbig); 
	//
	//
	// <x,x>
	TMatrixDSub(BigCov, 0, 2, 0, 2) = fXvCov; // <x.x>
	// 
	for(Int_t i=0; i<fNtr; i++){	// loop on <x, p_i>
		// <x,p>
		TMatrixD CovPX = BuildPXcov(i);
		TMatrixD CovXP(TMatrixD::kTransposed, CovPX);
		TMatrixDSub(BigCov, 0, 2, 3 * (i + 1), 3 * (i + 2) - 1) = CovXP; // <x.p>
		TMatrixDSub(BigCov, 3 * (i + 1), 3 * (i + 2) - 1, 0, 2) = CovPX; // <p.x>
		// <p, p>
		for(Int_t j=0; j<fNtr; j++){
			TMatrixD CovPP = BuildPPcov(i,j);
			TMatrixDSub(BigCov, 3 * (i + 1), 3 * (i + 2) - 1, 3 * (j + 1), 3 * (j + 2) - 1) = CovPP;
		}
	}
	// Symmetrize
	for(Int_t i=0; i<Nbig; i++){
		for(Int_t j=0; j<Nbig; j++)fBigCov(i,j) = 0.5*(BigCov(i,j)+BigCov(j,i));
	}
}
//
void VertexMore::FillBigPar()
{
	fBigPar.SetSub(0,fXv);
	for(Int_t i=0; i<fNtr; i++){
		TVector3 p = GetMomentum(i);
		Double_t p_arr[3] = {p.x(), p.y(), p.z()};
		TVectorD pv(3,p_arr);
		fBigPar.SetSub(3*i+3,pv);  
	}
}
//
// Vertex parameter errors
//
TMatrixD VertexMore::DparDx(TVector3 xv, TVector3 pv, Double_t Q)
{
	//
	// Derivative of track parameters wrt point on  track
	//
	// Track parameters
	TVectorD Par(5);
	if(Q != 0.0){				// Charged
		if(fUnits) xv *= 1.0e-3;	// Change to meters
		Par = XPtoPar(xv, pv, Q);
		if(fUnits) xv *= 1.0e3;
		if(fUnits)Par = ParToMm(Par);	// Back to mm
	}
	else Par = XPtoPar_N(xv, pv);		// Neutral
		
	//
	Double_t D = Par(0);
	Double_t ph0 = Par(1);
	Double_t lm = Par(4);
	//
	// Derivative matrix
	TMatrixD dParX(5,3); dParX.Zero();
	if(Q != 0.0){
		Double_t a = fa * Q;
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
		dParX(3, 0) = pv.Z()*dParX(1, 0)/a;
		dParX(3, 1) = pv.Z()*dParX(1, 1)/a;
		dParX(3, 2) = 1.0;
	}
	else{
		Double_t pt = Par(2);
		// D derivatives
		Double_t cs = TMath::Cos(ph0);
		Double_t sn = TMath::Sin(ph0);
		Double_t s = xv.Y()*sn+xv.X()*cs;
		//std::cout<<"dParX calc: sn= "<<sn<<", pt= "<<pt<<", D= "<<D<<", s= "<<s<<std::endl;
		dParX(0,0) = -sn;	// x
		dParX(0,1) =  cs;	// y
		// z0 derivatives
		dParX(3,0) = -lm*cs;	// x
		dParX(3,1) = -lm*sn;	// y
		dParX(3,2) = 1.0;	// z
	}
	//
	return dParX;
}
//
TMatrixD VertexMore::DparDp(TVector3 xv, TVector3 pv, Double_t Q)
{
	//
	// Derivative of track parameters wrt track momentum
	//
	// Track parameters
	TVectorD Par(5);
	if(Q != 0.0){				// Charged
		if(fUnits) xv *= 1.0e-3;	// Change to meters
		Par = XPtoPar(xv, pv, Q);
		if(fUnits) xv *= 1.0e3;
		if(fUnits)Par = ParToMm(Par);	// Back to mm
	}
	else Par = XPtoPar_N(xv, pv);		// Neutral	
	//
	Double_t D = Par(0);
	Double_t ph0 = Par(1);
	Double_t lm = Par(4);
	//
	// Derivative matrix
	TMatrixD dParP(5, 3); dParP.Zero();
	if(Q != 0.0){
		Double_t a = fa * Q;
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
		Double_t dsdpx = -pv.Y()/(pt*pt)-dParP(1,0);
		Double_t dsdpy =  pv.X()/(pt*pt)-dParP(1,1); 
		Double_t s = TMath::ATan2(pv.Y(),pv.X()) - ph0;
		if(s >  TMath::Pi()) s-= TMath::TwoPi();
		if(s < -TMath::Pi()) s+= TMath::TwoPi();
		dParP(3, 0) = -pv.Z()*dsdpx/a;
		dParP(3, 1) = -pv.Z()*dsdpy/a;
		dParP(3, 2) = -s/a; 
	}
	else{
		Double_t cs = TMath::Cos(ph0);
		Double_t sn = TMath::Sin(ph0);
		Double_t pt = pv.Pt();
		Double_t s = xv.Y()*sn+xv.X()*cs;
		//std::cout<<"dParP calc: sn= "<<sn<<", pt= "<<pt<<", D= "<<D<<", s= "<<s<<std::endl;
		// D derivatives
		dParP(0, 0) =  s*sn/pt;	// px
		dParP(0, 1) = -s*cs/pt;	// py
		// Phi0 derivatives
		dParP(1, 0) = -sn/pt;	// px
		dParP(1, 1) =  cs/pt;	// py
		// pt derivatives
		dParP(2, 0) = cs;	// px
		dParP(2, 1) = sn;	// py
		// z0 derivatives
		dParP(3, 0) = lm*(s*cs+D*sn)/pt;	// px
		dParP(3, 1) = lm*(s*sn-D*cs)/pt;	// py
		dParP(3, 2) = -s/pt;			// pz
		// ctg derivatives
		dParP(4, 0) = -lm*cs/pt;		// px
		dParP(4, 1) = -lm*sn/pt;		// py
		dParP(4, 2) = 1.0/pt;			// pz 
	}
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
	//
	// Vertex position
	TVector3 xv(fXv(0), fXv(1), fXv(2));
	TVector3 pv = fP;
	//
	// Parameter derivatives dPar/dx, dPar/dp
	TMatrixD dParX = DparDx(xv, pv, fQtot);	// (5, 3)
	TMatrixD dParP = DparDp(xv, pv, fQtot); // (5, 3)
	//
	TMatrixD dAlfdq(5,3*(fNtr+1));
	fCov.Zero();
	//
	TMatrixDSub(dAlfdq, 0, 4, 0, 2) = dParX;
	for(Int_t i=0; i<fNtr; i++)TMatrixDSub(dAlfdq, 0, 4, 3*(i+1), 3*(i+2)-1) = dParP;
	TMatrixDSym Bm = fBigCov;
	//if (!CheckPosDef(Bm)) std::cout<<"Error making BigCov"<<std::endl;
	fCov = Bm.Similarity(dAlfdq);
	if(!CheckPosDef(fCov))std::cout<<"VertexMore:: Error making fCov"<<std::endl;
	//	
	//
	return fCov;
}
//
//	Mass constraints
//
void VertexMore::AddMassConstraint(Double_t Mass, Int_t Nt, Double_t* masses, Int_t* list)
{	
	//
	fNc++;					// Update nr. of constraints
	fMassC.ResizeTo(fNc);			// Constraint vector
	fMderC.ResizeTo(fNc,3*(fNtr+1));	// Constraint derivatives wrt momenta
	//
	fMass.push_back(Mass);			// Mass of ith constraint
	fNtracks.push_back(Nt);			// nr. track in ith constraint
	fmasses.push_back(masses);		// List of masses of track list
	flists.push_back(list);			// Input list of tracks
	//
}

//
// Updated mass constraints and their derivatives
void VertexMore::UpdateConstr(TVectorD pp)
{
	//
	if(fNc < 1) std::cout<<"No mass constraints available. Do nothing."<<std::endl;
	else{
		fMassC.Zero();
		fMderC.Zero();
		for(Int_t n=0; n<fNc; n++){
			//
			// Invariant masses of constraints
			TLorentzVector pmu(0., 0., 0., 0.);	// Total momentum of nth constraint
			for(Int_t i=0; i<fNtracks[n]; i++){	// Loop on constraint tracks
				Int_t k = flists[n][i];
				Double_t mi = fmasses[n][i];
				TVectorD pin = pp.GetSub(3*k+3,3*k+5);
				TVector3 pi(pin(0),pin(1), pin(2));
				Double_t Ei = TMath::Sqrt(mi*mi+pi.Mag2());
				TLorentzVector pk(pi,Ei);		// 4-momentum of kth track 
				pmu += pk;
			}
			Double_t M = pmu.M();
			Double_t E = pmu.E();
			TVector3 p = pmu.Vect();
			fMassC[n] = M-fMass[n];			// Invariant mass
			//
			// Derivatives of invariant masses
			TVectorD pdfull(3*(fNtr+1)); pdfull.Zero();
			for(Int_t i=0; i<fNtracks[n]; i++){	// Loop on constraint tracks
				Int_t k = flists[n][i];
				Double_t mi = fmasses[n][i];
				TVectorD pin = pp.GetSub(3*k+3,3*k+5);
				TVector3 pi(pin(0), pin(1), pin(2));
				Double_t Ei = TMath::Sqrt(mi*mi+pi.Mag2());
				TVector3 pder = (1./M)*((E/Ei)*pi-p);
				Double_t pdrv[3]; pder.GetXYZ(pdrv);TVectorD pd(3,pdrv);
				pdfull.SetSub(3*k+3, pd);
			}
			TMatrixDRow(fMderC,n) = pdfull;
		}
	}
}
//
// Trigger mass constraint fit
void VertexMore::MassConstrFit()
{
	TVectorD p0 = fBigPar;		// Starting values
	TVectorD pp = p0;		// Initial expansion point
	TMatrixDSym S = fBigCov;	// Error matrix;
	TVectorD p = fBigPar; p.Zero();
	TMatrixDSym CovP = fBigCov; CovP.Zero();
	//
	Double_t eps = 1.0e-9;		// Fit precision
	Double_t deps = 1000.;		// Starting accuracy
	Int_t Nmax = 100;		// Maximum iterations
	Int_t Niter = 0;
	//
	while(TMath::Abs(deps) > eps){	// Main fit iteration loop
		UpdateConstr(pp);	// Update constraint and derivatives
		TMatrixD B = fMderC;	// Derivative matrix
		TMatrixD Bt(TMatrixD::kTransposed,fMderC);
		TVectorD f0 = fMassC;	// Constraint vector
		TMatrixDSym S0 = S;
		TMatrixDSym Wm1 = S0.Similarity(B);
		TMatrixDSym W = RegInv(Wm1);
		TMatrixD Ui(TMatrixD::kUnit, S);
		TMatrixDSym BtWB = W;
		BtWB.Similarity(Bt);
		TVectorD dp0 = p0-pp;
		TVectorD dp = (Ui-S*BtWB)*dp0;
		dp -=    (S*(Bt*W))*f0;
		TMatrixDSym SBtWBS = BtWB;
		SBtWBS.Similarity(S);
		CovP = S-SBtWBS;
		//
		TMatrixDSym Cinv = RegInv(S);
		deps = Cinv.Similarity(dp);
		//
		p = dp+pp;
		pp = p;
		Niter++;
		//
		if(Niter>Nmax) std::cout<<"VertexMore::MassConstrFit Maximum iterations reached"<<std::endl;
		if(Niter>Nmax) break;
	}
	//
	// Update everything
	//
	fBigPar = p;
	fBigCov = CovP;
	// Vertex
	fXv = p.GetSub(0,2);
	fXvCov = CovP.GetSub(0, 2, 0, 2);
	// Momenta
	TVector3 Ptot(0.,0.,0);
	for(Int_t i=0; i<fNtr; i++){
		TVectorD pm = p.GetSub(3*i+3,3*i+5);
		TVector3 p3(pm(0), pm(1), pm(2));
		Ptot += p3;
                if ( fpi[i] ) delete fpi[i];
		fpi[i] = new TVector3(p3);
		TMatrixDSym Cpm = CovP.GetSub(3*i+3,3*i+5,3*i+3,3*i+5);
                if ( fCpi[i] ) delete fCpi[i] ;
		fCpi[i] = new TMatrixDSym(Cpm);
	}
	fP = Ptot;
	TMatrixD PtotCov(3,3); PtotCov.Zero();
	for(Int_t i=0; i<fNtr; i++){
		for(Int_t k=0; k<fNtr; k++){
			TMatrixD Block = CovP.GetSub(3*i+3,3*i+5,3*k+3,3*k+5);
			PtotCov += Block;
		}
	}
	for(Int_t i=0; i<3; i++){
		for(Int_t k=0; k<3; k++)fCp(i,k) = 0.5*(PtotCov(i,k)+PtotCov(k,i));
	}
	//
	// Vertex track
	fPar = MakeVpar();
	fCov = MakeVcov();
	
}
