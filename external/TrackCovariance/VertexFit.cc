#include <TMath.h>
#include <TVectorD.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include "VertexFit.h"
//
// Constructors
//
//
// Empty construction (to be used when adding tracks later with AddTrk() )
VertexFit::VertexFit()
{
	fNtr = 0;
	fRold = -1.0;
	fVtxDone = kFALSE;
	fVtxCst = kFALSE;
	fxCst.ResizeTo(3);
	fCovCst.ResizeTo(3, 3);
	fCovCstInv.ResizeTo(3, 3);
	fXv.ResizeTo(3);
	fcovXv.ResizeTo(3, 3);
}
//
// Build from list of parameters and covariances
VertexFit::VertexFit(Int_t Ntr, TVectorD** trkPar, TMatrixDSym** trkCov)
{
	fNtr = Ntr;
	fRold = -1.0;
	fVtxDone = kFALSE;
	fVtxCst = kFALSE;
	fxCst.ResizeTo(3);
	fCovCst.ResizeTo(3, 3);
	fCovCstInv.ResizeTo(3, 3);
	fXv.ResizeTo(3);
	fcovXv.ResizeTo(3, 3);
	//
	for (Int_t i = 0; i < fNtr; i++)
	{
		TVectorD pr = *trkPar[i];
		fPar.push_back(new TVectorD(pr));
		fParNew.push_back(new TVectorD(pr));
		TMatrixDSym cv = *trkCov[i];
		fCov.push_back(new TMatrixDSym(cv));
		fCovNew.push_back(new TMatrixDSym(cv));
	}
	fChi2List.ResizeTo(fNtr);
	//
}
//
// Build from ObsTrk list of tracks
VertexFit::VertexFit(Int_t Ntr, ObsTrk** track)
{
	fNtr = Ntr;
	fRold = -1.0;
	fVtxDone = kFALSE;
	fVtxCst = kFALSE;
	fxCst.ResizeTo(3);
	fCovCst.ResizeTo(3, 3);
	fCovCstInv.ResizeTo(3, 3);
	fXv.ResizeTo(3);
	fcovXv.ResizeTo(3, 3);
	//
	fChi2List.ResizeTo(fNtr);
	for (Int_t i = 0; i < fNtr; i++)
	{
		fPar.push_back(new TVectorD(track[i]->GetObsPar()));
		fParNew.push_back(new TVectorD(track[i]->GetObsPar()));
		fCov.push_back(new TMatrixDSym(track[i]->GetCov()));
		fCovNew.push_back(new TMatrixDSym(track[i]->GetCov()));
	}
}
//
// Destructor
//
void VertexFit::ResetWrkArrays()
{
	Int_t N = (Int_t)ffi.size();
	for (Int_t i = 0; i < N; i++)
	{
		if (fx0i[i])  { fx0i[i]->Clear();		delete fx0i[i]; }
		if (fai[i])   { fai[i]->Clear();		delete fai[i]; }
		if (fdi[i])   { fdi[i]->Clear();		delete fdi[i]; }
		if (fAti[i])  { fAti[i]->Clear();		delete fAti[i]; }
		if (fDi[i])   { fDi[i]->Clear();		delete fDi[i]; }
		if (fWi[i])   { fWi[i]->Clear();		delete fWi[i]; }
		if (fWinvi[i]){ fWinvi[i]->Clear();     delete fWinvi[i]; }
	}
	fa2i.clear();
	fx0i.clear();
	fai.clear();
	fdi.clear();
	fAti.clear();
	fDi.clear();
	fWi.clear();
	fWinvi.clear();
}
VertexFit::~VertexFit()
{
	fxCst.Clear();		
	fCovCst.Clear();
	fCovCstInv.Clear();
	fXv.Clear();		
	fcovXv.Clear();		
	fChi2List.Clear();	
	//
	for (Int_t i = 0; i < fNtr; i++)
	{
		fPar[i]->Clear();		delete fPar[i];
		fParNew[i]->Clear();	delete fParNew[i];
		fCov[i]->Clear();		delete fCov[i];
		fCovNew[i]->Clear();	delete fCovNew[i];
	}
	fPar.clear();
	fParNew.clear();
	fCov.clear();
	fCovNew.clear();		
	//
	ResetWrkArrays();
	ffi.clear();
	fNtr = 0;
}
//
Double_t VertexFit::FastRv(TVectorD p1, TVectorD p2)
{
	//
	// Find radius of minimum distance between two tracks
	//
	// p = (D,phi, C, z0, ct)
	//
	// Define arrays
	//
	Double_t C1 = p1(2);
	Double_t C2 = p2(2);
	Double_t ph1 = p1(1);
	Double_t ph2 = p2(1);
	TVectorD x0 = Fill_x0(p1);
	TVectorD y0 = Fill_x0(p2);
	TVectorD n = Fill_a(p1, 0.0);
	n *= (2*C1);
	TVectorD k = Fill_a(p2, 0.0);
	k *= (2*C2);
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
	//
	// Higher order corrections
	X(0) += -C1 * smin(0) * smin(0) * TMath::Sin(ph1);
	X(1) +=  C1 * smin(0) * smin(0) * TMath::Cos(ph1);
	Y(0) += -C2 * smin(1) * smin(1) * TMath::Sin(ph2);
	Y(1) +=  C2 * smin(1) * smin(1) * TMath::Cos(ph2);
	//
	TVectorD Xavg = 0.5 * (X + Y);
	//
	//
	return TMath::Sqrt(Xavg(0)*Xavg(0)+Xavg(1)*Xavg(1));
}
//
// Starting radius determination
Double_t VertexFit::StartRadius()
{
	//
	// Maximum impact parameter
	Double_t Rd = 0;
	for (Int_t i = 0; i < fNtr; i++)
	{
		TVectorD par = *fPar[i];
		Double_t Dabs = TMath::Abs(par(0));
		if (Dabs > Rd)Rd = Dabs;
	}
	//-----------------------------
	//
	// Find track pair with phi difference closest to pi/2
	Int_t isel = 0; Int_t jsel = 0;		// selected track indices
	Double_t dSinMax = 0.0;				// Max phi difference
	for (Int_t i = 0; i < fNtr - 1; i++)
	{
		TVectorD pari = *fPar[i];
		Double_t phi1 = pari(1);

		for (Int_t j = i + 1; j < fNtr; j++)
		{
			TVectorD parj = *fPar[j];
			Double_t phi2 = parj(1);
			Double_t Sindphi = TMath::Abs(TMath::Sin(phi2 - phi1));
			if (Sindphi > dSinMax)
			{
				isel = i; jsel = j;
				dSinMax = Sindphi;
			}
		}
	}
	//
	//------------------------------------------
	//
	// Find radius of minimum distrance between tracks
	TVectorD p1 = *fPar[isel];
	TVectorD p2 = *fPar[jsel];
	Double_t R = FastRv(p1, p2);
	//
	R = 0.9 * R + 0.1 * Rd;		// Protect for overshoot
	//
	return R;
}
//
// Regularized symmetric matrix inversion
//
TMatrixDSym VertexFit::RegInv(TMatrixDSym& Min)
{
	TMatrixDSym M = Min;				// Decouple from input
	Int_t N = M.GetNrows();			// Matrix size
	TMatrixDSym D(N); D.Zero();		// Normaliztion matrix
	TMatrixDSym R(N);				// Normarized matrix
	TMatrixDSym Rinv(N);				// Inverse of R
	TMatrixDSym Minv(N);				// Inverse of M
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
void VertexFit::UpdateTrkArrays(Int_t i)
{
	//
	// Get track parameters and their covariance
	TVectorD par = *fParNew[i];
	TMatrixDSym Cov = *fCov[i];
	//
	// Fill all track related work arrays arrays
	Double_t fs = ffi[i];						// Get phase
	TVectorD xs = Fill_x(par, fs);
	fx0i.push_back(new TVectorD(xs));			// Start helix position
	//
	TMatrixD A = Fill_A(par, fs);				// A = dx/da = derivatives wrt track parameters
	TMatrixD At(TMatrixD::kTransposed, A);		// A transposed
	fAti.push_back(new TMatrixD(At));			// Store A' 
	TVectorD di = A * (par - *fPar[i]);		// x-shift
	fdi.push_back(new TVectorD(di));			// Store x-shift
	TMatrixDSym Winv = Cov.Similarity(A);		// W^-1 = A*C*A'
	fWinvi.push_back(new TMatrixDSym(Winv));	// Store W^-1 matrix
	//
	TMatrixDSym W = RegInv(Winv);				// W = (A*C*A')^-1
	fWi.push_back(new TMatrixDSym(W));			// Store W matrix
	//
	TVectorD a = Fill_a(par, fs);				// a = dx/ds = derivatives wrt phase
	fai.push_back(new TVectorD(a));				// Store a
	//
	Double_t a2 = W.Similarity(a);
	fa2i.push_back(a2);							// Store a2
	//
	// Build D matrix
	TMatrixDSym B(3);
	B.Rank1Update(a, -1. / a2);
	B.Similarity(W);
	TMatrixDSym Ds = W + B;						// D matrix
	fDi.push_back(new TMatrixDSym(Ds));			// Store D matrix
}
//
void  VertexFit::VertexFitter()
{
	//std::cout << "VertexFitter: just in" << std::endl;
	if (fNtr < 2 && !fVtxCst){
		std::cout << "VertexFit::VertexFitter - Method called with less than 2 tracks - Aborting " << std::endl;
		std::exit(1);
	}
	//
	// Vertex fit
	//
	// Initial variable definitions
	TVectorD x(3);
	TMatrixDSym covX(3);
	Double_t Chi2 = 0;
	TVectorD x0 = fXv;	// If previous fit done
	if (fRold < 0.0) {
		// External constraint
		if (fVtxCst) fRold = TMath::Sqrt(fxCst(0) * fxCst(0) + fxCst(1) * fxCst(1));
		// No constraint
		else for (Int_t i = 0; i < 3; i++)x0(i) = 1000.;	// Set to arbitrary large value 
	}
	//
	// Starting vertex radius approximation
	//
	Double_t R = fRold;						// Use previous fit if available
	if (R < 0.0) R = StartRadius();			// Rough vertex estimate
	//std::cout << "Start radius: " << R << std::endl;
	//
	// Iteration properties
	//
	Int_t Ntry = 0;
	Int_t TryMax = 100;
	Double_t eps = 1.0e-9; // vertex stability
	Double_t epsi = 1000.;
	//
	// Iteration loop
	while (epsi > eps && Ntry < TryMax)		// Iterate until found vertex is stable
	{
		// Initialize arrays
		x.Zero();
		TVectorD cterm(3); TMatrixDSym H(3); TMatrixDSym DW1D(3);
		covX.Zero();		// Reset vertex covariance
		cterm.Zero();	// Reset constant term
		H.Zero();		// Reset H matrix
		DW1D.Zero();
		//
		// Reset work arrays
		//
		ResetWrkArrays();
		//
		// Start loop on tracks
		//
		for (Int_t i = 0; i < fNtr; i++)
		{
			// Get track helix parameters and their covariance matrix
			TVectorD par = *fPar[i];
			TMatrixDSym Cov = *fCov[i];
			//
			// For first iteration only
			Double_t fs;
			if (Ntry <= 0)	// Initialize all phases on first pass
			{
				Double_t D = par(0);
				Double_t C = par(2);
				Double_t arg = TMath::Max(1.0e-6, (R * R - D * D) / (1 + 2 * C * D));
				fs = 2 * TMath::ASin(C * TMath::Sqrt(arg));
				ffi.push_back(fs);
			}
			//
			// Update track related arrays
			//
			UpdateTrkArrays(i);
			TMatrixDSym Ds = *fDi[i];
			TMatrixDSym Winv = *fWinvi[i];
			TMatrixDSym DsW1Ds = Winv.Similarity(Ds);	// Service matrix to calculate covX
			//
			// Update global arrays
			DW1D += DsW1Ds;
			// Update hessian
			H += Ds;
			// update constant term
			TVectorD xs = *fx0i[i] - *fdi[i];
			TVectorD xx0 = *fx0i[i];
			/*
			std::cout << "Iter. " << Ntry << ", trk " << i << ", x= "
				<< xx0(0) << ", " << xx0(1) << ", " << xx0(2)<<
				", ph0= "<<par(1)<< std::endl;
			*/
			cterm += Ds * xs;
		}				// End loop on tracks
		// Some additions in case of external constraints
		if (fVtxCst) {
			H += fCovCstInv;
			cterm += fCovCstInv * fxCst;
			DW1D += fCovCstInv;
		}
		//
		// update vertex position
		TMatrixDSym H1 = RegInv(H);
		x = H1 * cterm;
		//
		// Update vertex covariance
		covX = DW1D.Similarity(H1);
		//
		// Update phases and chi^2
		Chi2 = 0.0;
		for (Int_t i = 0; i < fNtr; i++)
		{
			TVectorD lambda = (*fDi[i]) * (*fx0i[i] - x - *fdi[i]);
			TMatrixDSym Wm1 = *fWinvi[i];
			fChi2List(i) = Wm1.Similarity(lambda);
			Chi2 += fChi2List(i);
			TVectorD a = *fai[i];
			TVectorD b = (*fWi[i]) * (x - (*fx0i[i]));
			for (Int_t j = 0; j < 3; j++)ffi[i] += a(j) * b(j) / fa2i[i];
			TVectorD newPar = *fPar[i] - ((*fCov[i]) * (*fAti[i])) * lambda;
			fParNew[i] = new TVectorD(newPar);
		}
		// Add external constraint to Chi2
		if (fVtxCst) Chi2 += fCovCstInv.Similarity(x - fxCst);
		//
		TVectorD dx = x - x0;
		x0 = x;
		// update vertex stability
		TMatrixDSym Hess = RegInv(covX);
		epsi = Hess.Similarity(dx);
		Ntry++;
		//
		// Store result
		//
		fXv = x;			// Vertex position
		fcovXv = covX;		// Vertex covariance
		fChi2 = Chi2;		// Vertex fit Chi2
		//std::cout << "Found vertex " << fXv(0) << ", " << fXv(1) << ", " << fXv(2) << std::endl;
	}		// end of iteration loop
	//
	fVtxDone = kTRUE;		// Set fit completion flag
	fRold = TMath::Sqrt(fXv(0)*fXv(0) + fXv(1)*fXv(1));	// Store fit radius
	//
}
//
// Return fit vertex
TVectorD VertexFit::GetVtx()
{
	if (!fVtxDone) VertexFitter();
	return fXv;
}
//
// Return fit vertex covariance
TMatrixDSym VertexFit::GetVtxCov()
{
	if (!fVtxDone) VertexFitter();
	return fcovXv;
}
//
// Return fit vertex chi2
Double_t VertexFit::GetVtxChi2()
{
	if (!fVtxDone) VertexFitter();
	return fChi2;
}
//
// Return array of chi2 contributions from each track
TVectorD VertexFit::GetVtxChi2List()
{
	if (!fVtxDone) VertexFitter();
	return fChi2List;
}
//
// Handle tracks/constraints
void VertexFit::AddVtxConstraint(TVectorD xv, TMatrixDSym cov)	// Add gaussian vertex constraint
{
	//std::cout << "VertexFit::AddVtxConstraint: Not implemented yet" << std::endl;
	fVtxCst = kTRUE;				// Vertex constraint flag
	fxCst = xv;						// Constraint value
	fCovCst = cov;					// Constraint covariance
	fCovCstInv = cov;
	fCovCstInv.Invert();				// Its inverse
	//
	// Set starting vertex as external constraint
	fXv = fxCst;
	fcovXv = fCovCst;
}
//
// Adding tracks one by one
void VertexFit::AddTrk(TVectorD *par, TMatrixDSym *Cov)			// Add track to input list
{
	fNtr++;
	fChi2List.ResizeTo(fNtr);	// Resize chi2 array
	fPar.push_back(par);			// add new track
	fCov.push_back(Cov);
	fParNew.push_back(par);			// add new track
	fCovNew.push_back(Cov);
	//
	// Reset previous vertex temp arrays
	ResetWrkArrays();
	ffi.clear();
	fVtxDone = kFALSE;			// Reset vertex done flag
}
//
// Removing tracks one by one
void VertexFit::RemoveTrk(Int_t iTrk)	// Remove iTrk track
{
	fNtr--;
	fChi2List.Clear();
	fChi2List.ResizeTo(fNtr);		// Resize chi2 array
	fPar.erase(fPar.begin() + iTrk);		// Remove track
	fCov.erase(fCov.begin() + iTrk);
	fParNew.erase(fParNew.begin() + iTrk);		// Remove track
	fCovNew.erase(fCovNew.begin() + iTrk);
	//
	// Reset previous vertex temp arrays
	ResetWrkArrays();
	ffi.clear();
	fVtxDone = kFALSE;			// Reset vertex done flag
}
