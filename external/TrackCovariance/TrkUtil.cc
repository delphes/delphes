#include "TrkUtil.h"
#include <iostream>
#include <algorithm>
#include <TSpline.h>
#include <TDecompChol.h>

// Constructor
TrkUtil::TrkUtil(Double_t Bz)
{
	fBz = Bz;
	fGasSel = 0;				// Default is He-Isobuthane (90-10)
	fRmin = 0.0;				// Lower		DCH radius
	fRmax = 0.0;				// Higher	DCH radius
	fZmin = 0.0;				// Lower		DCH z
	fZmax = 0.0;				// Higher	DCH z
}
TrkUtil::TrkUtil()
{
	fBz = 2.0;				// Default is 2 Tesla
	fGasSel = 0;				// Default is He-Isobuthane (90-10)
	fRmin = 0.0;				// Lower		DCH radius
	fRmax = 0.0;				// Higher	DCH radius
	fZmin = 0.0;				// Lower		DCH z
	fZmax = 0.0;				// Higher	DCH z
}
//
// Destructor
TrkUtil::~TrkUtil()
{
	fBz = 0.0;
	fGasSel = 0;				// Default is He-Isobuthane (90-10)
	fRmin = 0.0;				// Lower		DCH radius
	fRmax = 0.0;				// Higher	DCH radius
	fZmin = 0.0;				// Lower		DCH z
	fZmax = 0.0;				// Higher	DCH z
}
//
// Distance between two lines
//
void TrkUtil::LineDistance(TVector3 x0, TVector3 y0, TVector3 dirx, TVector3 diry, Double_t &sx, Double_t &sy, Double_t &distance)
{
	TMatrixDSym M(2);
	M(0,0) = dirx.Mag2();
	M(1,1) = diry.Mag2();
	M(0,1) = -dirx.Dot(diry);
	M(1,0) = M(0,1);
	M.Invert();
	TVectorD c(2);
	c(0) = dirx.Dot(y0-x0);
	c(1) = diry.Dot(x0-y0);
	TVectorD st = M*c;
	//
	// Fill output
	sx = st(0);
	sy = st(1);
	//
	TVector3 x = x0+sx*dirx;
	TVector3 y = y0+sy*diry;
	TVector3 d = x-y;
	distance = d.Mag();
}
//
// Covariance smearing
//
TVectorD TrkUtil::CovSmear(TVectorD x, TMatrixDSym C)
{
	//
	// Check arrays
	//
	// Consistency of dimensions
	Int_t Nvec = x.GetNrows();
	Int_t Nmat = C.GetNrows();
	if (Nvec != Nmat || Nvec == 0)
	{
		std::cout << "TrkUtil::CovSmear: vector/matrix mismatch. Aborting." << std::endl;
		exit(EXIT_FAILURE);
	}
	// Positive diagonal elements
	for (Int_t i = 0; i < Nvec; i++)
	{
		if (C(i, i) <= 0.0)
		{
			std::cout << "TrkUtil::CovSmear: covariance matrix has negative diagonal elements. Aborting." << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	//
	// Do a Choleski decomposition and random number extraction, with appropriate stabilization
	//
	TMatrixDSym CvN = C;
	TMatrixDSym DCv(Nvec); DCv.Zero();
	TMatrixDSym DCvInv(Nvec); DCvInv.Zero();
	for (Int_t id = 0; id < Nvec; id++)
	{
		Double_t dVal = TMath::Sqrt(C(id, id));
		DCv(id, id) = dVal;
		DCvInv(id, id) = 1.0 / dVal;
	}
	CvN.Similarity(DCvInv);			// Normalize diagonal to 1
	TDecompChol Chl(CvN);
	Bool_t OK = Chl.Decompose();		// Choleski decomposition of normalized matrix
	if (!OK)
	{
		std::cout << "TrkUtil::CovSmear: covariance matrix is not positive definite. Aborting." << std::endl;
		exit(EXIT_FAILURE);
	}
	TMatrixD U = Chl.GetU();			// Get Upper triangular matrix
	TMatrixD Ut(TMatrixD::kTransposed, U); // Transposed of U (lower triangular)
	TVectorD r(Nvec);
	for (Int_t i = 0; i < Nvec; i++)r(i) = gRandom->Gaus(0.0, 1.0);		// Array of normal random numbers
	TVectorD xOut = x + DCv * (Ut * r);	// Observed parameter vector
	//
	return xOut;
}
//
// Helix parameters from position and momentum
// static
TVectorD TrkUtil::XPtoPar(TVector3 x, TVector3 p, Double_t Q, Double_t Bz)
{
	//
	TVectorD Par(5);
	// Transverse parameters
	Double_t a = -Q * Bz * cSpeed();			// Units are Tesla, GeV and meters
	Double_t pt = p.Pt();
	Double_t C = a / (2 * pt);			// Half curvature
	//std::cout << "ObsTrk::XPtoPar: fB = " << fB << ", a = " << a << ", pt = " << pt << ", C = " << C << std::endl;
	Double_t r2 = x(0) * x(0) + x(1) * x(1);
	Double_t cross = x(0) * p(1) - x(1) * p(0);
	Double_t T = TMath::Sqrt(pt * pt - 2 * a * cross + a * a * r2);
	Double_t phi0 = TMath::ATan2((p(1) - a * x(0)), (p(0) + a * x(1)));	// Phi0
	Double_t D;							// Impact parameter D
	if (pt < 10.0) D = (T - pt) / a;
	else D = (-2 * cross + a * r2) / (T + pt);
	//
	Par(0) = D;		// Store D
	Par(1) = phi0;	// Store phi0
	Par(2) = C;		// Store C
	//Longitudinal parameters
	Double_t ct = p(2) / pt;
	// Old
	/*
	Double_t B = C * TMath::Sqrt(TMath::Max(r2 - D * D, 0.0) / (1 + 2 * C * D));
	Double_t st = TMath::ASin(B) / C;
	Double_t z0;
	Double_t dot = x(0) * p(0) + x(1) * p(1);
	if (dot > 0.0) z0 = x(2) - ct * st;
	else z0 = x(2) + ct * st;
	*/
	// New
	Double_t s = TMath::ATan2(p(1),p(0)) - phi0;
	if(s >  TMath::Pi()) s-= TMath::TwoPi();
	if(s < -TMath::Pi()) s+= TMath::TwoPi();
	Double_t z0 = x(2)-ct*s/(2.*C);
	//
	Par(3) = z0;		// Store z0
	Par(4) = ct;		// Store cot(theta)
	//
	return Par;
}
// non-static
TVectorD TrkUtil::XPtoPar(TVector3 x, TVector3 p, Double_t Q)
{
	//
	TVectorD Par(5);
	Double_t Bz = fBz;
	Par = XPtoPar(x, p, Q, Bz);
	//
	return Par;
}
//
TVector3 TrkUtil::ParToX(TVectorD Par)
{
	Double_t D = Par(0);
	Double_t phi0 = Par(1);
	Double_t z0 = Par(3);
	//
	TVector3 Xval;
	Xval(0) = -D * sin(phi0);
	Xval(1) = D * cos(phi0);
	Xval(2) = z0;
	//
	return Xval;
}
//
TVector3 TrkUtil::ParToP(TVectorD Par)
{
	if (fBz == 0.0)std::cout << "TrkUtil::ParToP: Warning Bz not set" << std::endl;
	//
	return ParToP(Par, fBz);
}
//
TVector3 TrkUtil::ParToP(TVectorD Par, Double_t Bz)
{
	Double_t C = Par(2);
	Double_t phi0 = Par(1);
	Double_t ct = Par(4);
	//
	TVector3 Pval;
	Double_t pt = Bz * cSpeed() / TMath::Abs(2 * C);
	Pval(0) = pt * cos(phi0);
	Pval(1) = pt * sin(phi0);
	Pval(2) = pt * ct;
	//
	return Pval;
}
//
// Neutrals
//
//static
TVectorD TrkUtil::XPtoPar_N(TVector3 x, TVector3 p)
{
//
// Output neutral track parameter vector:
// (D, phi0, pt, z0, cot(theta))
	TVectorD pout(5);
// Pt
	pout(2) = p.Pt();
// Direction
	Double_t csp0 = p.X()/p.Pt();
	Double_t snp0 = p.Y()/p.Pt();
	pout(4) = p.Z()/p.Pt();
	pout(1) = TMath::ATan2(snp0,csp0);
// Impact parameters
	pout(0) = x.Y()*csp0-x.X()*snp0;	// D (transverse)
	Double_t s = x.Y()*snp0+x.X()*csp0;	// dist from pma
	pout(3) = x.Z()-pout(4)*s;		// Z0
	
//
	return pout;
}
//
// static
TVector3 TrkUtil::ParToP_N(TVectorD Par)
{
	Double_t phi0 = Par(1);
	Double_t pt = Par(2);
	Double_t ctg = Par(4);
//
	TVector3 p(pt*TMath::Cos(phi0), pt*TMath::Sin(phi0), pt*ctg);
	return p;
}
//
Double_t TrkUtil::ParToQ(TVectorD Par)
{
	return TMath::Sign(1.0, -Par(2));
}

//
// Parameter conversion to ACTS format
TVectorD TrkUtil::ParToACTS(TVectorD Par)
{
	TVectorD pACTS(6);	// Return vector
	//
	Double_t b = -cSpeed() * fBz / 2.;
	pACTS(0) = 1000 * Par(0);		// D from m to mm
	pACTS(1) = 1000 * Par(3);	// z0 from m to mm
	pACTS(2) = Par(1);			// Phi0 is unchanged
	pACTS(3) = atan2(1.0, Par(4));		// Theta in [0, pi] range
	pACTS(4) = Par(2) / (b * sqrt(1 + Par(4) * Par(4)));		// q/p in GeV
	pACTS(5) = 0.0;				// Time: currently undefined
	//
	return pACTS;
}
// Covariance conversion to ACTS format
TMatrixDSym TrkUtil::CovToACTS(TVectorD Par, TMatrixDSym Cov)
{
	TMatrixDSym cACTS(6); cACTS.Zero();
	Double_t b = -cSpeed() * fBz / 2.;
	//
	// Fill derivative matrix
	TMatrixD A(5, 5);	A.Zero();
	Double_t ct = Par(4);	// cot(theta)
	Double_t C = Par(2);		// half curvature
	A(0, 0) = 1000.;		// D-D	conversion to mm
	A(1, 2) = 1.0;		// phi0-phi0
	A(2, 4) = 1.0 / (sqrt(1.0 + ct * ct) * b);	// q/p-C
	A(3, 1) = 1000.;		// z0-z0 conversion to mm
	A(4, 3) = -1.0 / (1.0 + ct * ct); // theta - cot(theta)
	A(4, 4) = -C * ct / (b * pow(1.0 + ct * ct, 3.0 / 2.0)); // q/p-cot(theta)
	//
	TMatrixDSym Cv = Cov;
	TMatrixD At(5, 5);
	At.Transpose(A);
	Cv.Similarity(At);
	TMatrixDSub(cACTS, 0, 4, 0, 4) = Cv;
	cACTS(5, 5) = 0.1;	// Currently undefined: set to arbitrary value to avoid crashes
	//
	return cACTS;
}
//
// Parameter conversion to ILC format
TVectorD TrkUtil::ParToILC(TVectorD Par)
{
	TVectorD pILC(5);	// Return vector
	//
	pILC(0) = Par(0) * 1.0e3;			// d0 in mm
	pILC(1) = Par(1);				// phi0 is unchanged
	pILC(2) = -2 * Par(2) * 1.0e-3;	// w in mm^-1
	pILC(3) = Par(3) * 1.0e3;			// z0 in mm
	pILC(4) = Par(4);				// tan(lambda) = cot(theta)
	//
	return pILC;
}
// Covariance conversion to ILC format
TMatrixDSym TrkUtil::CovToILC(TMatrixDSym Cov)
{
	TMatrixDSym cILC(5); cILC.Zero();
	//
	// Fill derivative matrix
	TMatrixD A(5, 5);	A.Zero();
	//
	A(0, 0) = 1.0e3;		// D-d0 in mm
	A(1, 1) = 1.0;			// phi0-phi0
	A(2, 2) = -2.0e-3;		// w-C
	A(3, 3) = 1.0e3;		// z0-z0 conversion to mm
	A(4, 4) = 1.0;			// tan(lambda) - cot(theta)
	//
	TMatrixDSym Cv = Cov;
	TMatrixD At(5, 5);
	At.Transpose(A);
	Cv.Similarity(At);
	cILC = Cv;
	//
	return cILC;
}
//
// Conversion from meters to mm
TVectorD TrkUtil::ParToMm(TVectorD Par)				// Parameter conversion
{
	TVectorD Pmm(5);					// Return vector
	//
	Pmm(0) = Par(0) * 1.0e3;			// d0 in mm
	Pmm(1) = Par(1);					// phi0 is unchanged
	Pmm(2) = Par(2) * 1.0e-3;			// C in mm^-1
	Pmm(3) = Par(3) * 1.0e3;			// z0 in mm
	Pmm(4) = Par(4);					// tan(lambda) = cot(theta) unchanged
	//
	return Pmm;
}
TMatrixDSym TrkUtil::CovToMm(TMatrixDSym Cov)		// Covariance conversion
{
	TMatrixDSym Cmm(5); Cmm.Zero();
	//
	// Fill derivative matrix
	TMatrixD A(5, 5);	A.Zero();
	//
	A(0, 0) = 1.0e3;		// D-d0 in mm
	A(1, 1) = 1.0;			// phi0-phi0
	A(2, 2) = 1.0e-3;		// C-C
	A(3, 3) = 1.0e3;		// z0-z0 conversion to mm
	A(4, 4) = 1.0;			// lambda - cot(theta)
	//
	TMatrixDSym Cv = Cov;
	TMatrixD At(5, 5);
	At.Transpose(A);
	Cv.Similarity(At);
	Cmm = Cv;
	//
	return Cmm;
}//
// Regularized symmetric matrix inversion
//
TMatrixDSym TrkUtil::RegInv(TMatrixDSym& Min)
{
	TMatrixDSym M = Min;				// Decouple from input
	Int_t N = M.GetNrows();			// Matrix size
	TMatrixDSym D(N); D.Zero();		// Normaliztion matrix
	TMatrixDSym R(N);				// Normarized matrix
	TMatrixDSym Rinv(N);				// Inverse of R
	TMatrixDSym Minv(N);				// Inverse of M
	//
	//*******************
	// Trivial case N = 1
	//*******************
	//
	if(N == 1){
		Minv(0,0) = 1.0;
		if(M(0,0) != 0.0) Minv(0,0) = 1.0/M(0,0);
		return Minv;
	}
	//
	// Check for 0's and normalize
	for (Int_t i = 0; i < N; i++)
	{
		if (M(i, i) != 0.0) D(i, i) = 1. / TMath::Sqrt(TMath::Abs(M(i, i)));
		else D(i, i) = 1.0;
	}
	R = M.Similarity(D);
	//
	// Recursive algorithms stops when N = 2
	//
	//****************
	// case N = 2  ***
	//****************
	if (N == 2)
	{
		Double_t det = R(0, 0) * R(1, 1) - R(0, 1) * R(1, 0);
		if (det == 0)
		{
			std::cout << "VertexFit::RegInv: null determinant for N = 2" << std::endl;
			Rinv.Zero();	// Return null matrix
		}
		else
		{
			// invert matrix 
			Rinv(0, 0) = R(1, 1);
			Rinv(0, 1) = -R(0, 1);
			Rinv(1, 0) = Rinv(0, 1);
			Rinv(1, 1) = R(0, 0);
			Rinv *= 1. / det;
		}
	}
	//****************
	// case N > 2  ***
	//****************
	else
	{
		// Break up matrix
		TMatrixDSym Q = R.GetSub(0, N - 2, 0, N - 2);	// Upper left 
		TVectorD p(N - 1);
		for (Int_t i = 0; i < N - 1; i++)p(i) = R(N - 1, i);
		Double_t q = R(N - 1, N - 1);
		//Invert pieces and re-assemble
		TMatrixDSym Ainv(N - 1);
		TMatrixDSym A(N - 1);
		if (TMath::Abs(q) > 1.0e-15)
		{
			// Case |q| > 0
			Ainv.Rank1Update(p, -1.0 / q);
			Ainv += Q;
			A = RegInv(Ainv);		// Recursive call
			TMatrixDSub(Rinv, 0, N - 2, 0, N - 2) = A;
			//
			TVectorD b = (-1.0 / q) * (A * p);
			for (Int_t i = 0; i < N - 1; i++)
			{
				Rinv(N - 1, i) = b(i);
				Rinv(i, N - 1) = b(i);
			}
			//
			Double_t pdotb = 0.;
			for (Int_t i = 0; i < N - 1; i++)pdotb += p(i) * b(i);
			Double_t c = (1.0 - pdotb) / q;
			Rinv(N - 1, N - 1) = c;
		}
		else
		{
			// case q = 0
			TMatrixDSym Qinv = RegInv(Q);		// Recursive call
			Double_t a = Qinv.Similarity(p);
			Double_t c = -1.0 / a;
			Rinv(N - 1, N - 1) = c;
			//
			TVectorD b = (1.0 / a) * (Qinv * p);
			for (Int_t i = 0; i < N - 1; i++)
			{
				Rinv(N - 1, i) = b(i);
				Rinv(i, N - 1) = b(i);
			}
			//
			A.Rank1Update(p, -1 / a);
			A += Q;
			A.Similarity(Qinv);
			TMatrixDSub(Rinv, 0, N - 2, 0, N - 2) = A;
		}
	}
	Minv = Rinv.Similarity(D);
	return Minv;
}

//
// Check potive definite matrix
//

Bool_t TrkUtil::CheckPosDef(TMatrixDSym Msym)
{
	Bool_t retVal = kTRUE;
	Int_t N = Msym.GetNrows();
	//std::cout<<"N = "<<N<<", Msym = "; Msym.Print();
	TMatrixDSym Nsym(N);
	TVectorD Diag(N);
	for(Int_t i=0; i< N; i++){
		if(Msym(i,i) <= 0){
			std::cout<<"CheckDefPos: found <= 0 on main diagonal M("<<i<<", "<<i<<") = "<<
			Msym(i,i)<<std::endl;
			retVal = kFALSE;
		}
		else Diag(i) = TMath::Sqrt(Msym(i,i));
	}
	//
	if(retVal){
	//
	// Normalize
		for(Int_t i=0; i<N; i++){
			for(Int_t j=0; j<N; j++)Nsym(i,j) = Msym(i,j)/(Diag(i)*Diag(j));
		}
	}
	//
	// Find eigenvalues
	//
	TMatrixDSymEigen Eign(Nsym);
	TVectorD lambda = Eign.GetEigenValues();
	for(Int_t i=0; i< N; i++){
		if(lambda(i) <= 0){
			std::cout<<"CheckDefPos: found <= 0 eigenvalue E("<<i<<") = "<<
			lambda(i)<<std::endl;
			std::cout<<"CheckDefPos: input matrix NOT posite definite. Printing normalized matrix."<<std::endl;
			Nsym.Print();
			retVal = kFALSE;
		}
	}
	//
	return retVal;
}

//
// Track tracjectory
//
TVector3 TrkUtil::Xtrack(TVectorD par, Double_t s)
{
	//
	// unpack parameters
	Double_t D = par(0);
	Double_t p0 = par(1);
	Double_t C = par(2);
	Double_t z0 = par(3);
	Double_t ct = par(4);
	//
	Double_t x = -D * TMath::Sin(p0) + (TMath::Sin(s + p0) - TMath::Sin(p0)) / (2 * C);
	Double_t y =  D * TMath::Cos(p0) - (TMath::Cos(s + p0) - TMath::Cos(p0)) / (2 * C);	
	Double_t z = z0 + ct * s / (2 * C);
	//
	TVector3 Xt(x, y, z);
	return Xt;
}
//
// Phase
//
Double_t TrkUtil::GetPhase(TVectorD x, TVectorD par)
{
	// Definitions
	// Transverse track parameters
	Double_t D = par(0);
	Double_t phi0 = par(1);
	Double_t sf = TMath::Sin(phi0);
	Double_t cf = TMath::Cos(phi0);
	Double_t C = par(2);
	//
	Double_t sins = 2.*C*(x(0)*cf+x(1)*sf);
	Double_t s = TMath::ASin(sins);
	//
	return s;
}
//
//	Phase derivatives
//	Track parameters
TVectorD TrkUtil::dsdPar(TVectorD x, TVectorD par)
{
	// 
	// Definitions
	// Transverse track parameters
	Double_t D = par(0);
	Double_t phi0 = par(1);
	Double_t sf = TMath::Sin(phi0);
	Double_t cf = TMath::Cos(phi0);
	Double_t C = par(2);
	//
	Double_t sins = 2.*C*(x(0)*cf+x(1)*sf);
	Double_t coss = TMath::Sqrt(1.-sins*sins);
	//
	// Derivatives
	//
	TVectorD der(5); der.Zero();
	der(1) = (2.*C*(-x(0)*sf+x(1)*cf))/coss;
	der(2) = 2.*(x(0)*cf+x(1)*sf)/coss;
//
	return der;
}
//
// position
TVectorD TrkUtil::dsdx(TVectorD x, TVectorD par)
{
	// 
	// Definitions
	// Transverse track parameters
	Double_t D = par(0);
	Double_t phi0 = par(1);
	Double_t sf = TMath::Sin(phi0);
	Double_t cf = TMath::Cos(phi0);
	Double_t C = par(2);
	//
	Double_t sins = 2.*C*(x(0)*cf+x(1)*sf);
	Double_t coss = TMath::Sqrt(1.-sins*sins);
	//
	// Derivatives
	//
	TVectorD der(3); der.Zero();
	der(0) = 2.*C*cf/coss;
	der(1) = 2.*C*sf/coss;
//
	return der;
}
//
// Trajectory of neutrals
//
TVector3 TrkUtil::Xtrack_N(TVectorD par, Double_t s)
{
	Double_t p0 = par(1);
	Double_t ctg = par(4);
	TVector3 x0 = ParToX(par);
	TVector3 dir(TMath::Cos(p0), TMath::Sin(p0), ctg);
	TVector3 Xt = x0 + s*dir;
//
	return Xt;
}
//
// Track derivatives
//
// Constant radius
// R-Phi
TVectorD TrkUtil::derRphi_R(TVectorD par, Double_t R)
{
	TVectorD dRphi(5);	// return vector
	//
	// unpack parameters
	Double_t D = par(0);
	Double_t C = par(2);
	//
	Double_t s = 2 * TMath::ASin(C * TMath::Sqrt((R * R - D * D)/(1 + 2 * C * D)));
	TVector3 X = Xtrack(par, s);		// Intersection point
	TVector3 v(-X.y()/R, X.x()/R, 0.);	// measurement direction
	TMatrixD derX = derXdPar(par, s);	// dX/dp
	TVectorD derXs = derXds(par, s);	// dX/ds
	TVectorD ders = dsdPar_R(par, R);	// ds/dp	
	//
	for (Int_t i = 0; i < 5; i++)
	{
		dRphi(i) = 0.;
		for (Int_t j = 0; j < 3; j++)
		{
			dRphi(i) += v(j) * (derX(j, i) + derXs(j) * ders(i));
		}
	}
	//
	return dRphi;
}
// z
TVectorD TrkUtil::derZ_R(TVectorD par, Double_t R)
{

	TVectorD dZ(5);	// return vector
	//
	// unpack parameters
	Double_t D = par(0);
	Double_t C = par(2);
	//
	Double_t s = 2 * TMath::ASin(C * TMath::Sqrt((R * R - D * D)/(1 + 2 * C * D))); // phase
	TVector3 v(0., 0., 1.);				// measurement direction
	TMatrixD derX = derXdPar(par, s);	// dX/dp
	TVectorD derXs = derXds(par, s);	// dX/ds
	TVectorD ders = dsdPar_R(par, R);	// ds/dp	
	//
	for (Int_t i = 0; i < 5; i++)
	{
		dZ(i) = 0.;
		for (Int_t j = 0; j < 3; j++)
		{
			dZ(i) += v(j) * (derX(j, i) + derXs(j) * ders(i));
		}
	}
	//
	return dZ;
}
//
// constant z
// R-Phi
TVectorD TrkUtil::derRphi_Z(TVectorD par, Double_t z)
{
	TVectorD dRphi(5);	// return vector
	//
	// unpack parameters
	Double_t C = par(2);
	Double_t z0 = par(3);
	Double_t ct = par(4);
	//
	Double_t s = 2 * C * (z - z0) / ct;
	TVector3 X = Xtrack(par, s);			// Intersection point
	TVector3 v(-X.y() / X.Pt(), X.x() / X.Pt(), 0.);	// measurement direction
	TMatrixD derX = derXdPar(par, s);		// dX/dp
	TVectorD derXs = derXds(par, s);		// dX/ds
	TVectorD ders = dsdPar_z(par, z);		// ds/dp	
	//
	for (Int_t i = 0; i < 5; i++)
	{
		dRphi(i) = 0.;
		for (Int_t j = 0; j < 3; j++)
		{
			dRphi(i) += v(j) * (derX(j, i) + derXs(j) * ders(i));
		}
	}
	//
	return dRphi;

}
// R
TVectorD TrkUtil::derR_Z(TVectorD par, Double_t z)
{
	TVectorD dR(5);	// return vector
	//
	// unpack parameters
	Double_t C = par(2);
	Double_t z0 = par(3);
	Double_t ct = par(4);
	//
	Double_t s = 2 * C * (z - z0) / ct;
	TVector3 X = Xtrack(par, s);			// Intersection point
	TVector3 v(X.x() / X.Pt(), X.y() / X.Pt(), 0.);	// measurement direction
	TMatrixD derX = derXdPar(par, s);		// dX/dp
	TVectorD derXs = derXds(par, s);		// dX/ds
	TVectorD ders = dsdPar_z(par, z);	// ds/dp	
	//
	for (Int_t i = 0; i < 5; i++)
	{
		dR(i) = 0.;
		for (Int_t j = 0; j < 3; j++)
		{
			dR(i) += v(j) * (derX(j, i) + derXs(j) * ders(i));
		}
	}
	//
	return dR;

}
//
// derivatives of track trajectory
//
// dX/dPar
TMatrixD TrkUtil::derXdPar(TVectorD par, Double_t s)
{
	TMatrixD dxdp(3, 5);	// return matrix
	//
	// unpack parameters
	Double_t D = par(0);
	Double_t p0 = par(1);
	Double_t C = par(2);
	//Double_t z0 = par(3);
	Double_t ct = par(4);
	//
	// derivatives
	// dx/dD
	dxdp(0, 0) = -TMath::Sin(p0);
	dxdp(1, 0) =  TMath::Cos(p0);
	dxdp(2, 0) = 0.;
	// dx/dphi0
	dxdp(0, 1) = -D * TMath::Cos(p0) + (TMath::Cos(s + p0) - TMath::Cos(p0)) / (2 * C);
	dxdp(1, 1) = -D * TMath::Sin(p0) + (TMath::Sin(s + p0) - TMath::Sin(p0)) / (2 * C);
	dxdp(2, 1) = 0;
	// dx/dC
	dxdp(0, 2) = -(TMath::Sin(s + p0) - TMath::Sin(p0)) / (2 * C * C);
	dxdp(1, 2) =  (TMath::Cos(s + p0) - TMath::Cos(p0)) / (2 * C * C);
	dxdp(2, 2) = -ct * s / (2 * C * C);
	// dx/dz0
	dxdp(0, 3) = 0;
	dxdp(1, 3) = 0;
	dxdp(2, 3) = 1.;
	// dx/dCtg
	dxdp(0, 4) = 0;
	dxdp(1, 4) = 0;
	dxdp(2, 4) = s / (2 * C);
	//
	return dxdp;
}
//
// dX/ds
//
TVectorD TrkUtil::derXds(TVectorD par, Double_t s)
{
	TVectorD dxds(3);	// return vector
	//
	// unpack parameters
	Double_t p0 = par(1);
	Double_t C = par(2);
	Double_t ct = par(4);
	//
	// dX/ds
	dxds(0) = TMath::Cos(s + p0) / (2 * C);
	dxds(1) = TMath::Sin(s + p0) / (2 * C);
	dxds(2) = ct / (2 * C);
	//
	return dxds;
}
//
// derivative of trajectory phase s
//Constant R
TVectorD TrkUtil::dsdPar_R(TVectorD par, Double_t R)
{
	TVectorD dsdp(5);	// return vector
	//
	// unpack parameters
	Double_t D = par(0);
	//Double_t p0 = par(1);
	Double_t C = par(2);
	//
	// derivatives
	Double_t opCD = 1. + 2 * C * D;
	Double_t A = C*TMath::Sqrt((R*R-D*D)/opCD);
	Double_t sqA0 = TMath::Sqrt(1. - A * A);
	Double_t dMin = 0.01;
	Double_t sqA = TMath::Max(dMin, sqA0);	// Protect against divergence
	//
	dsdp(0) = -2 * C * C * (D * (1. + C * D) + C * R * R) / (A * sqA * opCD * opCD);
	dsdp(1) = 0;
	dsdp(2) = 2 * A * (1 + C * D) / (C * sqA * opCD);
	dsdp(3) = 0;
	dsdp(4) = 0;
	//
	return dsdp;
}
// Constant z
TVectorD TrkUtil::dsdPar_z(TVectorD par, Double_t z)
{
	TVectorD dsdp(5);	// return vector
	//
	// unpack parameters
	Double_t C = par(2);
	Double_t z0 = par(3);
	Double_t ct = par(4);
	//
	// derivatives
	//
	dsdp(0) = 0;
	dsdp(1) = 0;
	dsdp(2) = 2*(z-z0)/ct;
	dsdp(3) = -2*C/ct;
	dsdp(4) = -2*C*(z-z0)/(ct*ct);
	//
	return dsdp;
}
//
// Derivatives of neutral trajectory
//dX/dPar
TMatrixD TrkUtil::derXdPar_N(TVectorD par, Double_t s)	// derivatives of position wrt parameters
{
	TMatrixD dxdp(3, 5);	// return matrix
	//
	// unpack parameters
	Double_t D = par(0);
	Double_t p0 = par(1);
	//Double_t pt = par(2);
	//Double_t z0 = par(3);
	//Double_t ct = par(4);
	//
	//
	// derivatives
	// dx/dD
	dxdp(0, 0) = -TMath::Sin(p0);
	dxdp(1, 0) =  TMath::Cos(p0);
	dxdp(2, 0) = 0.;
	// dx/dphi0
	dxdp(0, 1) = -D * TMath::Cos(p0) - s * TMath::Sin(p0);
	dxdp(1, 1) = -D * TMath::Sin(p0) + s * TMath::Cos(p0);
	dxdp(2, 1) = 0;
	// dx/dpt
	dxdp(0, 2) = 0.;
	dxdp(1, 2) = 0.;
	dxdp(2, 2) = 0.;
	// dx/dz0
	dxdp(0, 3) = 0;
	dxdp(1, 3) = 0;
	dxdp(2, 3) = 1.;
	// dx/dCtg
	dxdp(0, 4) = 0;
	dxdp(1, 4) = 0;
	dxdp(2, 4) = s;
//
	return dxdp;
}
//dX/ds 
TVectorD TrkUtil::derXds_N(TVectorD par, Double_t s)	// derivatives of position wrt phase
{
	TVectorD dxds(3);	// return vector
	//
	// unpack parameters
	Double_t p0 = par(1);
	Double_t ct = par(4);
	//
	// dX/ds
	dxds(0) = TMath::Cos(p0);
	dxds(1) = TMath::Sin(p0);
	dxds(2) = ct;
	//
	return dxds;
}
//ds/dPar const R
TVectorD TrkUtil::dsdPar_R_N(TVectorD par, Double_t R)	// derivatives of phase at constant R
{
	TVectorD dsdp(5);	// return vector
	//
	// unpack parameters
	Double_t D = par(0);
	//
	// derivatives
	Double_t s = TMath::Sqrt(R*R-D*D);
	Double_t dMin = 0.01;
	Double_t deriv = TMath::Max(dMin, -D/s);	// Protect against divergence
	//
	dsdp(0) = deriv;
	dsdp(1) = 0;
	dsdp(2) = 0;
	dsdp(3) = 0;
	dsdp(4) = 0;
	//
	return dsdp;
}
//ds/dPar const z
TVectorD TrkUtil::dsdPar_z_N(TVectorD par, Double_t z)	// derivatives of phase at constant z
{
	TVectorD dsdp(5);	// return vector
	//
	// unpack parameters
	Double_t z0 = par(3);
	Double_t ct = par(4);
	//
	// derivatives
	//
	dsdp(0) = 0;
	dsdp(1) = 0;
	dsdp(2) = 0;
	dsdp(3) = -1./ct;
	dsdp(4) = -(z-z0)/(ct*ct);
	//
	return dsdp;
}
//
// Setup chamber volume
void TrkUtil::SetDchBoundaries(Double_t Rmin, Double_t Rmax, Double_t Zmin, Double_t Zmax)
{
	fRmin = Rmin;				// Lower		DCH radius
	fRmax = Rmax;				// Higher	DCH radius
	fZmin = Zmin;				// Lower		DCH z
	fZmax = Zmax;				// Higher	DCH z
}
//
// Get Trakck length inside DCH volume
Double_t TrkUtil::TrkLen(TVectorD Par) const
{
	Double_t tLength = 0.0;
	// Check if geometry is initialized
	if (fZmin == 0.0 && fZmax == 0.0)
	{
		// No geometry set so send a warning and return 0
		std::cout << "TrkUtil::TrkLen() called without a DCH volume defined" << std::endl;
	}
	else
	{
		//******************************************************************
		// Determine the track length inside the chamber   ****
		//******************************************************************
		//
		// Track pararameters
		Double_t D = Par(0);		// Transverse impact parameter
		//Double_t phi0 = Par(1);		// Transverse direction at minimum approach
		Double_t C = Par(2);		// Half curvature
		Double_t z0 = Par(3);		// Z at minimum approach
		Double_t ct = Par(4);		// cot(theta)
		//std::cout << "TrkUtil:: parameters: D= " << D << ", phi0= " << phi0
		//	<< ", C= " << C << ", z0= " << z0 << ", ct= " << ct << std::endl;
		//
		// Track length per unit phase change 
		Double_t Scale = sqrt(1.0 + ct * ct) / (2.0 * TMath::Abs(C));
		//
		// Find intersections with chamber boundaries
		//
		Double_t phRin = 0.0;			// phase of inner cylinder 
		Double_t phRin2 = 0.0;			// phase of inner cylinder intersection (2nd branch)
		Double_t phRhi = 0.0;			// phase of outer cylinder intersection
		Double_t phZmn = 0.0;			// phase of left wall intersection
		Double_t phZmx = 0.0;			// phase of right wall intersection
		//  ... with inner cylinder
		Double_t Rtop = TMath::Abs((1.0 + C * D) / C);

		if (Rtop > fRmin && TMath::Abs(D) < fRmin) // *** don't treat large D tracks for the moment ***
		{
			Double_t ph = 2 * asin(C * sqrt((fRmin * fRmin - D * D) / (1.0 + 2.0 * C * D)));
			Double_t z = z0 + ct * ph / (2.0 * C);

			//std::cout << "Rin intersection: ph = " << ph<<", z= "<<z << std::endl;

			if (z < fZmax && z > fZmin)	phRin = TMath::Abs(ph);	// Intersection inside chamber volume	
			//
			// Include second branch of loopers
			Double_t Pi = 3.14159265358979323846;
			Double_t ph2 = 2 * Pi - TMath::Abs(ph);
			if (ph < 0)ph2 = -ph2;
			z = z0 + ct * ph2 / (2.0 * C);
			if (z < fZmax && z > fZmin)	phRin2 = TMath::Abs(ph2);	// Intersection inside chamber volume
		}
		//  ... with outer cylinder
		if (Rtop > fRmax && TMath::Abs(D) < fRmax) // *** don't treat large D tracks for the moment ***
		{
			Double_t ph = 2 * asin(C * sqrt((fRmax * fRmax - D * D) / (1.0 + 2.0 * C * D)));
			Double_t z = z0 + ct * ph / (2.0 * C);
			if (z < fZmax && z > fZmin)	phRhi = TMath::Abs(ph);	// Intersection inside chamber volume	
		}
		//  ... with left wall
		Double_t Zdir = (fZmin - z0) / ct;
		if (Zdir > 0.0)
		{
			Double_t ph = 2.0 * C * Zdir;
			Double_t Rint = sqrt(D * D + (1.0 + 2.0 * C * D) * pow(sin(ph / 2), 2) / (C * C));
			if (Rint < fRmax && Rint > fRmin)	phZmn = TMath::Abs(ph);	// Intersection inside chamber volume	
		}
		//  ... with right wall
		Zdir = (fZmax - z0) / ct;
		if (Zdir > 0.0)
		{
			Double_t ph = 2.0 * C * Zdir;
			Double_t Rint = sqrt(D * D + (1.0 + 2.0 * C * D) * pow(sin(ph / 2), 2) / (C * C));
			if (Rint < fRmax && Rint > fRmin)	phZmx = TMath::Abs(ph);	// Intersection inside chamber volume	
		}
		//
		// Order phases and keep the lowest two non-zero ones
		//
		const Int_t Nint = 5;
		Double_t dPhase = 0.0;	// Phase difference between two close intersections
		Double_t ph_arr[Nint] = { phRin, phRin2, phRhi, phZmn, phZmx };
		std::sort(ph_arr, ph_arr + Nint);
		Int_t iPos = -1;		// First element > 0
		for (Int_t i = 0; i < Nint; i++)
		{
			if (ph_arr[i] <= 0.0) iPos = i;
		}

		if (iPos < Nint - 2)
		{
			dPhase = ph_arr[iPos + 2] - ph_arr[iPos + 1];
			tLength = dPhase * Scale;
		}
	}
	return tLength;
}
//
// Return number of ionization clusters
Bool_t TrkUtil::IonClusters(Double_t& Ncl, Double_t mass, TVectorD Par)
{
	//
	// Units are meters/Tesla/GeV
	//
	Ncl = 0.0;
	Bool_t Signal = kFALSE;
	Double_t tLen = 0;
	// Check if geometry is initialized
	if (fZmin == 0.0 && fZmax == 0.0)
	{
		// No geometry set so send a warning and return 0
		std::cout << "TrkUtil::IonClusters() called without a volume defined" << std::endl;
	}
	else tLen = TrkLen(Par);

	//******************************************************************
	// Now get the number of clusters                       ****
	//******************************************************************
	//
	Double_t muClu = 0.0;	// mean number of clusters
	Double_t bg = 0.0;		// beta*gamma
	Ncl = 0.0;
	if (tLen > 0.0)
	{
		Signal = kTRUE;
		//
		// Find beta*gamma
		if (fBz == 0.0)
		{
			Signal = kFALSE;
			std::cout << "TrkUtil::IonClusters: Please set Bz!!!" << std::endl;
		}
		else
		{
			TVector3 p = ParToP(Par);
			bg = p.Mag() / mass;
			muClu = Nclusters(bg) * tLen;				// Avg. number of clusters

			Ncl = gRandom->PoissonD(muClu);			// Actual number of clusters
		}

	}
	//
	return Signal;
}
//
//
Double_t TrkUtil::Nclusters(Double_t begam)
{
	Int_t Opt = fGasSel;
	Double_t Nclu = Nclusters(begam, Opt);
	//
	return Nclu;
}
//
Double_t TrkUtil::Nclusters(Double_t begam, Int_t Opt) {
	//
	// Opt = 0: He 90 - Isobutane 10
	//     = 1: pure He
	//     = 2: Argon 50 - Ethane 50
	//     = 3: pure Argon
	//
	//
	const Int_t Npt = 19;
	Double_t bg[Npt] = { 0.5, 0.8, 1., 2., 3., 4., 5., 8., 10.,
	12., 15., 20., 50., 100., 200., 500., 1000., 10000., 20000.};
	//
	// He 90 - Isobutane 10
	Double_t ncl_He_Iso[Npt] = { 42.94, 23.6,18.97,12.98,12.2,12.13,
	12.24,12.73,13.03,13.29,13.63,14.08,15.56,16.43,16.8,16.95,16.98, 16.98, 16.98};
	//
	// pure He
	Double_t ncl_He[Npt] = { 11.79,6.5,5.23,3.59,3.38,3.37,3.4,3.54,3.63,
				3.7,3.8,3.92,4.33,4.61,4.78,4.87,4.89, 4.89, 4.89};
	//
	// Argon 50 - Ethane 50
	Double_t ncl_Ar_Eth[Npt] = { 130.04,71.55,57.56,39.44,37.08,36.9,
	37.25,38.76,39.68,40.49,41.53,42.91,46.8,48.09,48.59,48.85,48.93,48.93, 48.93};
	//
	// pure Argon
	Double_t ncl_Ar[Npt] = { 88.69,48.93,39.41,27.09,25.51,25.43,25.69,
	26.78,27.44,28.02,28.77,29.78,32.67,33.75,34.24,34.57,34.68, 34.68, 34.68};
	//
	Double_t ncl[Npt];
	switch (Opt)
	{
	case 0: std::copy(ncl_He_Iso, ncl_He_Iso + Npt, ncl);	// He-Isobutane
		break;
	case 1: std::copy(ncl_He, ncl_He + Npt, ncl);		// pure He
		break;
	case 2: std::copy(ncl_Ar_Eth, ncl_Ar_Eth + Npt, ncl);	// Argon - Ethane
		break;
	case 3: std::copy(ncl_Ar, ncl_Ar + Npt, ncl);		// pure Argon
		break;
	}
	//
	Double_t interp = 0.0;
	TSpline3* sp3 = new TSpline3("sp3", bg, ncl, Npt);
	if (begam > bg[0] && begam < bg[Npt - 1]) interp = sp3->Eval(begam);
	return 100 * interp;
}
//
Double_t TrkUtil::funcNcl(Double_t* xp, Double_t* par) {
	Double_t bg = xp[0];
	return Nclusters(bg);
}
//
void TrkUtil::SetGasMix(Int_t Opt)
{
	if (Opt < 0 || Opt > 3)
	{
		std::cout << "TrkUtil::SetGasMix Gas option not allowed. No action."
			<< std::endl;
	}
	else fGasSel = Opt;
}
