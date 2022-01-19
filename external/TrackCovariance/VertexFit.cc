/*
Vertex fitting code
*/
#include <TMath.h>
#include <TVectorD.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
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
		fPar.push_back(new TVectorD(*trkPar[i]));
		fParNew.push_back(new TVectorD(*trkPar[i]));
		fCov.push_back(new TMatrixDSym(*trkCov[i]));
		fCovNew.push_back(new TMatrixDSym(*trkCov[i]));
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
	Int_t N = (Int_t)fdi.size();
	if(N > 0){
		for (Int_t i = 0; i < N; i++)
		{
			if (fx0i[i])  { fx0i[i]->Clear();	delete fx0i[i]; }
			if (fai[i])   { fai[i]->Clear();	delete fai[i]; }
			if (fdi[i])   { fdi[i]->Clear();	delete fdi[i]; }
			if (fAti[i])  { fAti[i]->Clear();	delete fAti[i]; }
			if (fDi[i])   { fDi[i]->Clear();	delete fDi[i]; }
			if (fWi[i])   { fWi[i]->Clear();	delete fWi[i]; }
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
//
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
	// Get track parameters, covariance and phase
	Double_t fs = ffi[i];			// Get phase
	TVectorD par = *fParNew[i];
	TMatrixDSym Cov = *fCov[i];
	//
	// Fill all track related work arrays arrays
	TMatrixD A = derXdPar(par, fs);		// A = dx/da = derivatives wrt track parameters
	TMatrixDSym Winv = Cov;				
	Winv.Similarity(A);			// W^-1 = A*C*A'

	TMatrixD At(TMatrixD::kTransposed, A);		// A transposed
	fAti.push_back(new TMatrixD(At));		// Store A'
	fWinvi.push_back(new TMatrixDSym(Winv));	// Store W^-1 matrix
	//
	TVectorD xs = Fill_x(par, fs);
	fx0i.push_back(new TVectorD(xs));			// Start helix position
	// 
	TVectorD di = A * (par - *fPar[i]);			// x-shift
	fdi.push_back(new TVectorD(di));			// Store x-shift	
	//
	TMatrixDSym W = RegInv(Winv);				// W = (A*C*A')^-1
	fWi.push_back(new TMatrixDSym(W));			// Store W matrix
	//
	TVectorD a = derXds(par, fs);				// a = dx/ds = derivatives wrt phase
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
void VertexFit::VtxFitNoSteer()
{
	//
	// Initialize
	//
	std::vector<TVectorD*> x0i;				// Tracks at ma
	std::vector<TVectorD*> ni;				// Track derivative wrt phase
	std::vector<TMatrixDSym*> Ci;			// Position error matrix at fixed phase
	std::vector<TVectorD*> wi;				// Ci*ni
	//
	// Track loop
	for (Int_t i = 0; i < fNtr; i++)
	{
		Double_t s = 0.;
		TVectorD par = *fPar[i];
		TMatrixDSym Cov = *fCov[i];
		x0i.push_back(new TVectorD(Fill_x0(par)));
		ni.push_back(new TVectorD(derXds(par, s)));
		TMatrixD A = derXdPar(par, s);
		Ci.push_back(new TMatrixDSym(Cov.Similarity(A)));
		TMatrixDSym Cinv = RegInv(*Ci[i]);
		wi.push_back(new TVectorD(Cinv * (*ni[i])));
	}
	//std::cout << "Vtx init completed. fNtr = "<<fNtr << std::endl;
	//
	// Get fit vertex
	//
	TMatrixDSym D(3); D.Zero();
	TVectorD Dx(3); Dx.Zero();
	for (Int_t i = 0; i < fNtr; i++)
	{
		TMatrixDSym Cinv = RegInv(*Ci[i]);
		TMatrixDSym W(3);
		W.Rank1Update(*wi[i], 1. / Ci[i]->Similarity(*wi[i]));
		TMatrixDSym Dd = Cinv - W;
		D += Dd;
		Dx += Dd * (*x0i[i]);
	}
	if(fVtxCst){
		D  += fCovCstInv;
		Dx += fCovCstInv*fxCst;
	}
	fXv = RegInv(D) * Dx;
	//std::cout << "Fast vertex (x, y, z) = "<<fXv(0)<<", "<<fXv(1)<<", "<<fXv(2) << std::endl;
	//
	// Get fit phases
	//
	for (Int_t i = 0; i < fNtr; i++){
		Double_t si = Dot(*wi[i], fXv - (*x0i[i])) / Ci[i]->Similarity(*wi[i]);
		ffi.push_back(si);
		//TVectorD xvi = Fill_x(*fPar[i],si);
		//std::cout << "Fast vertex "<<i<<": xvi = "<<xvi(0)<<", "<<xvi(1)<<", "<<xvi(2) 
		//	<<", si = "<<si<< std::endl;
	}
	//
	// Cleanup
	for (Int_t i = 0; i < fNtr; i++)
	{
		x0i[i]->Clear();	delete x0i[i];
		ni[i]->Clear();		delete ni[i];
		Ci[i]->Clear();		delete Ci[i];
		wi[i]->Clear();		delete wi[i];
	}
	x0i.clear();;
	ni.clear();
	Ci.clear(); 
	wi.clear();
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
	//
	VtxFitNoSteer();	// Fast vertex finder on first pass (set ffi and fXv)
	TVectorD x0 = fXv;
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
		cterm.Zero();		// Reset constant term
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
			//TVectorD xx0 = *fx0i[i];
			
			//std::cout << "Iter. " << Ntry << ", trk " << i << ", xs= "
			//	<< xs(0) << ", " << xs(1) << ", " << xs(2)<<
			//	", ph0= "<<par(1)<< std::endl;
			
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
			if(fChi2List(i) < 0.0){
				//std::cout<<"# "<<i<<", Chi2= "<<fChi2List(i)<<", Wm1:"<<std::endl; Wm1.Print();
				//std::cout<<"Lambda= "<<std::endl; lambda.Print();
			}
			Chi2 += fChi2List(i);
			TVectorD a = *fai[i];
			TVectorD b = (*fWi[i]) * (x - *fx0i[i] + *fdi[i]);
			ffi[i] += Dot(a, b) / fa2i[i];
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
	fRold = TMath::Sqrt(fXv(0)*fXv(0) + fXv(1)*fXv(1));	// Store fit 
	//std::cout << "Found vertex " << fXv(0) << ", " << fXv(1) << ", " << fXv(2) 
	//	<< ", after "<<Ntry<<" iterations"<<std::endl;
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
