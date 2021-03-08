#include <TMath.h>
#include <TVectorD.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include "VertexFit.h"
//
// Constructors
//
// Empty
VertexFit::VertexFit()
{
	fNtr = 0;
	fVtxDone = kFALSE;
	fVtxCst = kFALSE;
	fxCst.ResizeTo(3);
	fCovCst.ResizeTo(3, 3);
	fXv.ResizeTo(3);
	fcovXv.ResizeTo(3, 3);
}
// Parameters and covariances
VertexFit::VertexFit(Int_t Ntr, TVectorD** trkPar, TMatrixDSym** trkCov)
{
	fNtr = Ntr;
	fVtxDone = kFALSE;
	fVtxCst = kFALSE;
	fxCst.ResizeTo(3);
	fCovCst.ResizeTo(3, 3);
	fXv.ResizeTo(3);
	fcovXv.ResizeTo(3, 3);
	//
	fPar = trkPar;
	fCov = trkCov;
	fChi2List.ResizeTo(Ntr);
	//
	ffi = new Double_t[Ntr];				// Fit phases
	fx0i = new TVectorD * [Ntr];			// Track expansion points
	for (Int_t i = 0; i < Ntr; i++) fx0i[i] = new TVectorD(3);
	fai = new TVectorD * [Ntr];			// dx/dphi
	for (Int_t i = 0; i < Ntr; i++) fai[i] = new TVectorD(3);
	fa2i = new Double_t[Ntr];				// a'Wa
	fDi = new TMatrixDSym * [Ntr];		// W-WBW
	for (Int_t i = 0; i < Ntr; i++) fDi[i] = new TMatrixDSym(3);
	fWi = new TMatrixDSym * [Ntr];	// (ACA')^-1
	for (Int_t i = 0; i < Ntr; i++) fWi[i] = new TMatrixDSym(3);
	fWinvi = new TMatrixDSym * [Ntr];	// ACA'
	for (Int_t i = 0; i < Ntr; i++) fWinvi[i] = new TMatrixDSym(3);
}
// ObsTrk list
VertexFit::VertexFit(Int_t Ntr, ObsTrk** track)
{
	fNtr = Ntr;
	fVtxDone = kFALSE;
	fVtxCst = kFALSE;
	fxCst.ResizeTo(3);
	fCovCst.ResizeTo(3, 3);
	fXv.ResizeTo(3);
	fcovXv.ResizeTo(3, 3);
	//
	fPar = new TVectorD * [Ntr];
	fCov = new TMatrixDSym * [Ntr];
	fChi2List.ResizeTo(Ntr);
	for (Int_t i = 0; i < Ntr; i++)
	{
		fPar[i] = new TVectorD(track[i]->GetObsPar());
		fCov[i] = new TMatrixDSym(track[i]->GetCov());
	}
	//
	ffi = new Double_t[Ntr];				// Fit phases
	fx0i = new TVectorD * [Ntr];			// Track expansion points
	for (Int_t i = 0; i < Ntr; i++) fx0i[i] = new TVectorD(3);
	fai = new TVectorD * [Ntr];			// dx/dphi
	for (Int_t i = 0; i < Ntr; i++) fai[i] = new TVectorD(3);
	fa2i = new Double_t[Ntr];				// a'Wa
	fDi = new TMatrixDSym * [Ntr];	// W-WBW
	for (Int_t i = 0; i < Ntr; i++) fDi[i] = new TMatrixDSym(3);
	fWi = new TMatrixDSym * [Ntr];	// (ACA')^-1
	for (Int_t i = 0; i < Ntr; i++) fWi[i] = new TMatrixDSym(3);
	fWinvi = new TMatrixDSym * [Ntr];	// ACA'
	for (Int_t i = 0; i < Ntr; i++) fWinvi[i] = new TMatrixDSym(3);
}
//
// Destructor
VertexFit::~VertexFit()
{
	fxCst.Clear();
	fCovCst.Clear();
	fXv.Clear();
	fcovXv.Clear();
	fChi2List.Clear();
	//
	for (Int_t i = 0; i < fNtr; i++)
	{
		fPar[i]->Clear();
		fCov[i]->Clear();
		//
		fx0i[i]->Clear(); delete fx0i[i];
		fai[i]->Clear(); delete fai[i];
		fDi[i]->Clear(); delete fDi[i];
		fWi[i]->Clear(); delete fWi[i];
		fWinvi[i]->Clear();	delete fWinvi[i];
	}
	fNtr = 0;
	delete[] fPar;
	delete[] fCov;
	delete[] ffi;
	delete[] fa2i;
	delete[] fx0i;
	delete[] fai;
	delete[] fDi;
	delete[] fWi;
	delete[] fWinvi;
}
//
Double_t VertexFit::FastRv1(TVectorD p1, TVectorD p2)
{
	//
	// Find radius of intersection between two tracks in the transverse plane
	//
	// p = (D,phi, C, z0, ct)
	//
	// Define arrays
	//
	Double_t r1 = 1.0 / p1(2);
	Double_t r2 = 1.0 / p2(2);
	TVectorD x0 = Fill_x0(p1);
	TVectorD y0 = Fill_x0(p2);
	TVectorD n = Fill_a(p1, 0.0);
	n *= r1;
	TVectorD k = Fill_a(p2, 0.0);
	k *= r2;
	//
	// Setup and solve linear system
	//
	Double_t nn = 0; for (Int_t i = 0; i < 3; i++)nn += n(i) * n(i);
	Double_t nk = 0; for (Int_t i = 0; i < 3; i++)nk += n(i) * k(i);
	Double_t kk = 0; for (Int_t i = 0; i < 3; i++)kk += k(i) * k(i);
	Double_t discr = nn * kk - nk * nk;
	TMatrixDSym H(2);
	H(0, 0) = kk;
	H(0, 1) = nk;
	H(1, 0) = H(0, 1);
	H(1, 1) = nn;
	TVectorD c(2);
	c(0) = 0; for (Int_t i = 0; i < 3; i++)c(0) += n(i) * (y0(i) - x0(i));
	c(1) = 0; for (Int_t i = 0; i < 3; i++)c(1) += -k(i) * (y0(i) - x0(i));
	TVectorD smin = (H * c);
	smin *= 1.0 / discr;
	//
	TVectorD X = x0 + smin(0) * n;
	TVectorD Y = y0 + smin(1) * k;
	Double_t R1 = TMath::Sqrt(X(0) * X(0) + X(1) * X(1));
	Double_t R2 = TMath::Sqrt(Y(0) * Y(0) + Y(1) * Y(1));
	//
	return 0.5 * (R1 + R2);
}
Double_t VertexFit::FastRv(TVectorD p1, TVectorD p2)
{
	//
	// Find radius of minimum distance
	//
	// p = (D,phi, C)
	//
	// Solving matrix
	TMatrixDSym H(2);
	H(0, 0) = -TMath::Cos(p2(1));
	H(0, 1) = TMath::Cos(p1(1));
	H(1, 0) = -TMath::Sin(p2(1));
	H(1, 1) = TMath::Sin(p1(1));
	Double_t Det = TMath::Sin(p2(1) - p1(1));
	H *= 1.0 / Det;
	//
	// Convergence parameters
	Int_t Ntry = 0;
	Int_t NtryMax = 100;
	Double_t eps = 1000.;
	Double_t epsMin = 1.0e-6;
	//
	// Vertex finding loop
	//
	TVectorD cterm(2);
	cterm(0) = p1(0);
	cterm(1) = p2(0);
	TVectorD xv(2);
	Double_t R = 1000.;
	while (eps > epsMin)
	{
		xv = H * cterm;
		Ntry++;
		if (Ntry > NtryMax)
		{
			std::cout << "FastRv: maximum number of iteration reached" << std::endl;
			break;
		}
		Double_t Rnew = TMath::Sqrt(xv(0) * xv(0) + xv(1) * xv(1));
		eps = Rnew - R;
		R = Rnew;
		cterm(0) = p1(2) * R * R;
		cterm(1) = p2(2) * R * R;
	}
	//
	return R;
}

TMatrixDSym VertexFit::RegInv3(TMatrixDSym& Smat0)
{
	//
	// Regularized inversion of symmetric 3x3 matrix with positive diagonal elements
	//
	TMatrixDSym Smat = Smat0;
	Int_t N = Smat.GetNrows();
	if (N != 3)
	{
		std::cout << "RegInv3 called with  matrix size != 3. Abort & return standard inversion." << std::endl;
		return Smat.Invert();
	}
	TMatrixDSym D(N); D.Zero();
	Bool_t dZero = kTRUE;	// No elements less or equal 0 on the diagonal
	for (Int_t i = 0; i < N; i++) if (Smat(i, i) <= 0.0)dZero = kFALSE;
	if (dZero)
	{
		for (Int_t i = 0; i < N; i++) D(i, i) = 1.0 / TMath::Sqrt(Smat(i, i));
		TMatrixDSym RegMat = Smat.Similarity(D);
		TMatrixDSym Q(2);
		for (Int_t i = 0; i < 2; i++)
		{
			for (Int_t j = 0; j < 2; j++)Q(i, j) = RegMat(i, j);
		}
		Double_t Det = 1 - Q(0, 1) * Q(1, 0);
		TMatrixDSym H(2);
		H = Q;
		H(0, 1) = -Q(0, 1);
		H(1, 0) = -Q(1, 0);
		TVectorD p(2);
		p(0) = RegMat(0, 2);
		p(1) = RegMat(1, 2);
		Double_t pHp = H.Similarity(p);
		Double_t h = pHp - Det;
		//
		TMatrixDSym pp(2); pp.Rank1Update(p);
		TMatrixDSym F = (h * H) - pp.Similarity(H);
		F *= 1.0 / Det;
		TVectorD b = H * p;
		TMatrixDSym InvReg(3);
		for (Int_t i = 0; i < 2; i++)
		{
			InvReg(i, 2) = b(i);
			InvReg(2, i) = b(i);
			for (Int_t j = 0; j < 2; j++) InvReg(i, j) = F(i, j);
		}
		InvReg(2, 2) = -Det;
		//
		InvReg *= 1.0 / h;
		//
		//
		return InvReg.Similarity(D);
	}
	else
	{
		D.Zero();
		for (Int_t i = 0; i < N; i++) D(i, i) = 1.0 / TMath::Sqrt(TMath::Abs(Smat(i, i)));
		TMatrixDSym RegMat = Smat.Similarity(D);
		RegMat.Invert();
		return RegMat.Similarity(D);
	}
}
//
//
//
TMatrixD VertexFit::Fill_A(TVectorD par, Double_t phi)
{
	//
	// Derivative of track 3D position vector with respect to track parameters at constant phase 
	//
	// par = vector of track parameters
	// phi = phase
	//
	TMatrixD A(3, 5);
	//
	// Decode input arrays
	//
	Double_t D = par(0);
	Double_t p0 = par(1);
	Double_t C = par(2);
	Double_t z0 = par(3);
	Double_t ct = par(4);
	//
	// Fill derivative matrix dx/d alpha
	// D
	A(0, 0) = -TMath::Sin(p0);
	A(1, 0) = TMath::Cos(p0);
	A(2, 0) = 0.0;
	// phi0
	A(0, 1) = -D * TMath::Cos(p0) + (TMath::Cos(phi + p0) - TMath::Cos(p0)) / (2 * C);
	A(1, 1) = -D * TMath::Sin(p0) + (TMath::Sin(phi + p0) - TMath::Sin(p0)) / (2 * C);
	A(2, 1) = 0.0;
	// C
	A(0, 2) = -(TMath::Sin(phi + p0) - TMath::Sin(p0)) / (2 * C * C);
	A(1, 2) = (TMath::Cos(phi + p0) - TMath::Cos(p0)) / (2 * C * C);
	A(2, 2) = -ct * phi / (2 * C * C);
	// z0
	A(0, 3) = 0.0;
	A(1, 3) = 0.0;
	A(2, 3) = 1.0;
	// ct = lambda
	A(0, 4) = 0.0;
	A(1, 4) = 0.0;
	A(2, 4) = phi / (2 * C);
	//
	return A;
}
//
TVectorD VertexFit::Fill_a(TVectorD par, Double_t phi)
{
	//
	// Derivative of track 3D position vector with respect to phase at constant track parameters
	//
	// par = vector of track parameters
	// phi = phase
	//
	TVectorD a(3);
	//
	// Decode input arrays
	//
	Double_t D = par(0);
	Double_t p0 = par(1);
	Double_t C = par(2);
	Double_t z0 = par(3);
	Double_t ct = par(4);
	//
	a(0) = TMath::Cos(phi + p0) / (2 * C);
	a(1) = TMath::Sin(phi + p0) / (2 * C);
	a(2) = ct / (2 * C);
	//
	return a;
}
//
TVectorD VertexFit::Fill_x0(TVectorD par)
{
	//
	// Calculate track 3D position at R = |D| (minimum approach to z-axis)
	//
	TVectorD x0(3);
	//
	// Decode input arrays
	//
	Double_t D = par(0);
	Double_t p0 = par(1);
	Double_t C = par(2);
	Double_t z0 = par(3);
	Double_t ct = par(4);
	//
	x0(0) = -D * TMath::Sin(p0);
	x0(1) = D * TMath::Cos(p0);
	x0(2) = z0;
	//
	return x0;
}
//
TVectorD VertexFit::Fill_x(TVectorD par, Double_t phi)
{
	//
	// Calculate track 3D position for a given phase, phi
	//
	TVectorD x(3);
	//
	// Decode input arrays
	//
	Double_t D = par(0);
	Double_t p0 = par(1);
	Double_t C = par(2);
	Double_t z0 = par(3);
	Double_t ct = par(4);
	//
	TVectorD x0 = Fill_x0(par);
	x(0) = x0(0) + (TMath::Sin(phi + p0) - TMath::Sin(p0)) / (2 * C);
	x(1) = x0(1) - (TMath::Cos(phi + p0) - TMath::Cos(p0)) / (2 * C);
	x(2) = x0(2) + ct * phi / (2 * C);
	//
	return x;
}
//
void  VertexFit::VertexFinder()
{
	//
	// Vertex fit (units are meters)
	//
	// Initial variable definitions
	TVectorD x(3);
	TMatrixDSym covX(3);
	TVectorD x0(3); for (Int_t v = 0; v < 3; v++)x0(v) = 100.; // set to large value
	Double_t Chi2 = 0;
	//
	// Stored quantities
	Double_t* fi = new Double_t[fNtr];				// Phases
	TVectorD** x0i = new TVectorD * [fNtr];			// Track expansion point
	TVectorD** ai = new TVectorD * [fNtr];				// dx/dphi
	Double_t* a2i = new Double_t[fNtr];				// a'Wa
	TMatrixDSym** Di = new TMatrixDSym * [fNtr];		// W-WBW
	TMatrixDSym** Wi = new TMatrixDSym * [fNtr];		// (ACA')^-1
	TMatrixDSym** Winvi = new TMatrixDSym * [fNtr];	// ACA'
	//
	// vertex radius approximation
	// Maximum impact parameter
	Double_t Rd = 0;
	for (Int_t i = 0; i < fNtr; i++)
	{
		//ObsTrk* t = tracks[i];
		TVectorD par = *fPar[i];
		Double_t Dabs = TMath::Abs(par(0));
		if (Dabs > Rd)Rd = Dabs;
	}
	//
	// Find track pair with largest phi difference
	Int_t isel = 0; Int_t jsel = 0; // selected track indices
	Double_t dphiMax = 0.0;	// Max phi difference

	for (Int_t i = 0; i < fNtr - 1; i++)
	{
		//ObsTrk* ti = tracks[i];
		TVectorD pari = *fPar[i];
		Double_t phi1 = pari(1);

		for (Int_t j = i + 1; j < fNtr; j++)
		{
			//ObsTrk* tj = tracks[j];
			TVectorD parj = *fPar[j];
			Double_t phi2 = parj(1);
			Double_t dphi = TMath::Abs(phi2 - phi1);
			if (dphi > TMath::Pi())dphi = TMath::TwoPi() - dphi;
			if (dphi > dphiMax)
			{
				isel = i; jsel = j;
				dphiMax = dphi;
			}
		}
	}
	//
	TVectorD p1 = *fPar[isel];
	TVectorD p2 = *fPar[jsel];
	Double_t R = FastRv1(p1, p2);
	if (R > 1000.0) R = Rd;
	R = 0.9 * R + 0.1 * Rd;
	//
	// Iteration properties
	//
	Int_t Ntry = 0;
	Int_t TryMax = 100;
	Double_t eps = 1.0e-9; // vertex stability
	Double_t epsi = 1000.;
	//
	while (epsi > eps && Ntry < TryMax)		// Iterate until found vertex is stable
	{
		x.Zero();
		TVectorD cterm(3); TMatrixDSym H(3); TMatrixDSym DW1D(3);
		covX.Zero();	// Reset vertex covariance
		cterm.Zero();	// Reset constant term
		H.Zero();		// Reset H matrix
		DW1D.Zero();
		//
		//std::cout << "VertexFinder: start loop on tracks" << std::endl;
		for (Int_t i = 0; i < fNtr; i++)
		{
			// Get track helix parameters and their covariance matrix
			TVectorD par = *fPar[i];
			TMatrixDSym Cov = *fCov[i];

			Double_t fs;
			if (Ntry <= 0)	// Initialize all phases on first pass
			{
				Double_t D = par(0);
				Double_t C = par(2);
				Double_t arg = TMath::Max(1.0e-6, (R * R - D * D) / (1 + 2 * C * D));
				fs = 2 * TMath::ASin(C * TMath::Sqrt(arg));
				fi[i] = fs;
			}
			//
			// Starting values
			//
			fs = fi[i];		// Get phase
			//std::cout << "VertexFinder: phase fs set" << std::endl;
			TVectorD xs = Fill_x(par, fs);
			//std::cout << "VertexFinder: position xs set" << std::endl;
			x0i[i] = new TVectorD(xs);				// Start helix position
			//std::cout << "VertexFinder: position x0i stored" << std::endl;
			// W matrix = (A*C*A')^-1; W^-1 = A*C*A'
			TMatrixD A = Fill_A(par, fs);			// A = dx/da = derivatives wrt track parameters
			//std::cout << "VertexFinder: derivatives A set" << std::endl;
			TMatrixDSym Winv = Cov.Similarity(A);	// W^-1 = A*C*A'
			Winvi[i] = new TMatrixDSym(Winv);		// Store W^-1 matrix
			//std::cout << "VertexFinder: Winvi stored" << std::endl;
			TMatrixDSym W = RegInv3(Winv);			// W = (A*C*A')^-1
			Wi[i] = new TMatrixDSym(W);				// Store W matrix
			//std::cout << "VertexFinder: Wi stored" << std::endl;
			TVectorD a = Fill_a(par, fs);			// a = dx/ds = derivatives wrt phase
			//std::cout << "VertexFinder: derivatives a set" << std::endl;
			ai[i] = new TVectorD(a);				// Store a
			//std::cout << "VertexFinder: derivatives a stored" << std::endl;
			Double_t a2 = W.Similarity(a);
			a2i[i] = a2;							// Store a2
			// Build D matrix
			TMatrixDSym B(3);

			B.Rank1Update(a, 1.0);
			B *= -1. / a2;
			B.Similarity(W);
			TMatrixDSym Ds = W + B;					// D matrix
			Di[i] = new TMatrixDSym(Ds);			// Store D matrix
			//std::cout << "VertexFinder: matrix Di stored" << std::endl;
			TMatrixDSym DsW1Ds = Winv.Similarity(Ds);	// Service matrix to calculate covX
			DW1D += DsW1Ds;
			// Update hessian
			H += Ds;
			// update constant term
			cterm += Ds * xs;
		}				// End loop on tracks
		//
		// update vertex position
		TMatrixDSym H1 = RegInv3(H);
		x = H1 * cterm;
		//std::cout << "VertexFinder: x vertex set" << std::endl;
		// Update vertex covariance
		covX = DW1D.Similarity(H1);
		//std::cout << "VertexFinder: cov vertex set" << std::endl;
		// Update phases and chi^2
		Chi2 = 0.0;
		for (Int_t i = 0; i < fNtr; i++)
		{
			TVectorD lambda = (*Di[i]) * (*x0i[i] - x);
			TMatrixDSym Wm1 = *Winvi[i];
			fChi2List(i) = Wm1.Similarity(lambda);
			Chi2 += fChi2List(i);
			TVectorD a = *ai[i];
			TVectorD b = (*Wi[i]) * (x - *x0i[i]);
			for (Int_t j = 0; j < 3; j++)fi[i] += a(j) * b(j) / a2i[i];
		}

		//
		TVectorD dx = x - x0;
		x0 = x;
		// update vertex stability
		TMatrixDSym Hess = RegInv3(covX);
		epsi = Hess.Similarity(dx);
		Ntry++;
		//
		// Store result
		//
		fXv = x;			// Vertex position
		fcovXv = covX;		// Vertex covariance
		fChi2 = Chi2;		// Vertex fit Chi2
		//
		// Store intermediate data
		//

		//std::cout << "VertexFinder: before store intermediate data" << std::endl;
		for (Int_t i = 0; i < fNtr; i++)
		{
			//std::cout << "VertexFinder: inside store intermediate data" << std::endl;
			//std::cout << "i = " << i << ", fi[i] = " << fi[i] << std::endl;
			//std::cout << "i = " << i << ", ffi[i] = " << ffi[i] << std::endl;
			ffi[i] = fi[i];				// Fit phases
			//std::cout << "VertexFinder: fi stored" << std::endl;
			fx0i[i] = x0i[i];			// Track expansion points
			//std::cout << "VertexFinder: x0i stored" << std::endl;
			fai[i] = ai[i];				// dx/dphi
			//std::cout << "VertexFinder: ai stored" << std::endl;
			fa2i[i] = a2i[i];			// a'Wa
			//std::cout << "VertexFinder: a2i stored" << std::endl;
			fDi[i] = Di[i];				// W-WBW
			//std::cout << "VertexFinder: Di stored" << std::endl;
			fWi[i] = Wi[i];				// (ACA')^-1
			//std::cout << "VertexFinder: Wi stored" << std::endl;
			fWinvi[i] = Winvi[i];		// ACA'
			//std::cout << "VertexFinder: Winvi stored" << std::endl;
		}
		//std::cout << "Iteration " << Ntry << " completed - Before cleanup" << std::endl;
		//
		// Cleanup
		//
		for (Int_t i = 0; i < fNtr; i++)
		{
			x0i[i]->Clear();
			Winvi[i]->Clear();
			Wi[i]->Clear();
			ai[i]->Clear();
			Di[i]->Clear();

			delete x0i[i];
			delete Winvi[i];
			delete Wi[i];
			delete ai[i];
			delete Di[i];
		}

		//std::cout << "Iteration " << Ntry << " completed - After cleanup" << std::endl;
	}
	//
	fVtxDone = kTRUE;	// Set fitting completion flag
	//
	delete[] fi;		// Phases
	delete[] x0i;		// Track expansion point
	delete[] ai;		// dx/dphi
	delete[] a2i;		// a'Wa
	delete[] Di;		// W-WBW
	delete[] Wi;		// (ACA')^-1
	delete[] Winvi;		// ACA'
}
//
TVectorD VertexFit::GetVtx()
{
	//std::cout << "GetVtx: flag set to " << fVtxDone << std::endl;
	if (!fVtxDone)VertexFinder();
	return fXv;
}

TMatrixDSym VertexFit::GetVtxCov()
{
	if (!fVtxDone)VertexFinder();
	return fcovXv;
}

Double_t VertexFit::GetVtxChi2()
{
	if (!fVtxDone)VertexFinder();
	return fChi2;
}

TVectorD VertexFit::GetVtxChi2List()
{
	if (!fVtxDone)VertexFinder();
	return fChi2List;
}
//
// Handle tracks/constraints
void VertexFit::AddVtxConstraint(TVectorD xv, TMatrixDSym cov)	// Add gaussian vertex constraint
{
	std::cout << "VertexFit::AddVtxConstraint: Not implemented yet" << std::endl;
}
//
void VertexFit::AddTrk(TVectorD par, TMatrixDSym Cov)			// Add track to input list
{
	std::cout << "VertexFit::AddTrk: Not implemented yet" << std::endl;
}
void VertexFit::RemoveTrk(Int_t iTrk)							// Remove iTrk track
{
	std::cout << "VertexFit::RemoveTrk: Not implemented yet" << std::endl;
}
