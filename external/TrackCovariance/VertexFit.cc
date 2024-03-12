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
	fRstart = -1.0;
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
	fRstart = -1.0;
	fVtxDone = kFALSE;
	fVtxCst = kFALSE;
	fxCst.ResizeTo(3);
	fCovCst.ResizeTo(3, 3);
	fCovCstInv.ResizeTo(3, 3);
	fXv.ResizeTo(3);
	fcovXv.ResizeTo(3, 3);
	//
	Bool_t Charged = kTRUE;
	for (Int_t i = 0; i < fNtr; i++)
	{
		fPar.push_back(new TVectorD(*trkPar[i]));
		fParNew.push_back(new TVectorD(*trkPar[i]));
		fCov.push_back(new TMatrixDSym(*trkCov[i]));
		fCovNew.push_back(new TMatrixDSym(*trkCov[i]));
		fCharged.push_back(Charged);
	}
	fChi2List.ResizeTo(fNtr);
	//
}
//
// Build from list of parameters and covariances with charge tag
VertexFit::VertexFit(Int_t Ntr, TVectorD** trkPar, TMatrixDSym** trkCov, Bool_t* Charged)
{
	fNtr = Ntr;
	fRstart = -1.0;
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
		fCharged.push_back(Charged[i]);
	}
	fChi2List.ResizeTo(fNtr);
	//
}
//
// Build from ObsTrk list of tracks
VertexFit::VertexFit(Int_t Ntr, ObsTrk** track)
{
	fNtr = Ntr;
	fRstart = -1.0;
	fVtxDone = kFALSE;
	fVtxCst = kFALSE;
	fxCst.ResizeTo(3);
	fCovCst.ResizeTo(3, 3);
	fCovCstInv.ResizeTo(3, 3);
	fXv.ResizeTo(3);
	fcovXv.ResizeTo(3, 3);
	//
	fChi2List.ResizeTo(fNtr);
	Bool_t Charged = kTRUE;
	for (Int_t i = 0; i < fNtr; i++)
	{
		fPar.push_back(new TVectorD(track[i]->GetObsPar()));
		fParNew.push_back(new TVectorD(track[i]->GetObsPar()));
		fCov.push_back(new TMatrixDSym(track[i]->GetCov()));
		fCovNew.push_back(new TMatrixDSym(track[i]->GetCov()));
		fCharged.push_back(Charged);
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
		fPar[i]->Clear();	delete fPar[i];
		fParNew[i]->Clear();	delete fParNew[i];
		fCov[i]->Clear();	delete fCov[i];
		fCovNew[i]->Clear();	delete fCovNew[i];
	}	
	fPar.clear();
	fParNew.clear();
	fCov.clear();
	fCovNew.clear();		
	//
	ResetWrkArrays();
	ffi.clear();	
	fCharged.clear();
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
	//Double_t C = par(2);
	Double_t z0 = par(3);
	//Double_t ct = par(4);
	//
	x0(0) = -D * TMath::Sin(p0);
	x0(1) = D * TMath::Cos(p0);
	x0(2) = z0;
	//
	return x0;
}
//
TVectorD VertexFit::Fill_x(TVectorD par, Double_t phi, Bool_t Charged)
{
	//
	// Calculate track 3D position for a given phase, phi
	//
	TVectorD x(3);
	TVector3 xt;
	if(Charged) xt = Xtrack(par, phi);
	else        xt = Xtrack_N(par,phi);
	for(Int_t i = 0; i<3; i++)x(i) = xt(i);
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
	TMatrixD A(3,5);				// A = dx/da = derivatives wrt track parameters
	if(fCharged[i]) A = derXdPar(par, fs);	
	else	        A = derXdPar_N(par, fs);
	TMatrixDSym Winv = Cov;				
	Winv.Similarity(A);			// W^-1 = A*C*A'

	TMatrixD At(TMatrixD::kTransposed, A);		// A transposed
	fAti.push_back(new TMatrixD(At));		// Store A'
	fWinvi.push_back(new TMatrixDSym(Winv));	// Store W^-1 matrix
	//
	TVectorD xs = Fill_x(par, fs, fCharged[i]);
	fx0i.push_back(new TVectorD(xs));			// Start helix position
	// 
	TVectorD di = A * (par - *fPar[i]);			// x-shift
	fdi.push_back(new TVectorD(di));			// Store x-shift	
	//
	TMatrixDSym W = RegInv(Winv);				// W = (A*C*A')^-1
	fWi.push_back(new TMatrixDSym(W));			// Store W matrix
	//
	TVectorD a(3);						// a = dx/ds = derivatives wrt phase
	if(fCharged[i]) a = derXds(par, fs);	
	else 		a = derXds_N(par, fs);	
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
	std::vector<TMatrixDSym*> Ci;				// Position error matrix at fixed phase
	std::vector<TVectorD*> wi;				// Ci*ni
	std::vector<Double_t> s_in;				// Starting phase
	//
	// 
	//
	// Track loop
	for (Int_t i = 0; i < fNtr; i++)
	{
		TVectorD par = *fPar[i];
		TMatrixDSym Cov = *fCov[i];
		Double_t s = 0.;
		// Case when starting radius is provided
		if(fRstart > TMath::Abs(par(0))){
			if(fCharged[i])s = 2.*TMath::ASin(par(2)*TMath::Sqrt((fRstart*fRstart-par(0)*par(0))/(1.+2.*par(2)*par(0))));
			else s = TMath::Sqrt(fRstart*fRstart-par(0)*par(0));
		}
		//
		x0i.push_back(new TVectorD(Fill_x(par, s, fCharged[i])));
		TMatrixD A(3,5);
		if(fCharged[i]){
			ni.push_back(new TVectorD(derXds(par, s)));
			A = derXdPar(par, s);}
		else{
			ni.push_back(new TVectorD(derXds_N(par, s)));
			A = derXdPar_N(par, s);}
		//
		Ci.push_back(new TMatrixDSym(Cov.Similarity(A)));
		TMatrixDSym Cinv = RegInv(*Ci[i]);
		wi.push_back(new TVectorD(Cinv * (*ni[i])));
		s_in.push_back(s);
		if(!fCharged[i]){
			/*std::cout<<"i = "<<i<<", Neutral s = "<<s<<std::endl;
			std::cout<<"Par = "; par.Print();
			std::cout<<"Cov = "; Cov.Print();
			std::cout<<"A: "; A.Print();
			std::cout<<"Ci "; Ci[i]->Print();
			std::cout<<"wi "; wi[i]->Print();
			*/
		}
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
		ffi.push_back(si+s_in[i]);
		//TVectorD xvi = Fill_x(*fPar[i],si, fCharged[i]);
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
	s_in.clear();
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
	Double_t eps = 1.0e-12; // vertex stability
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
                        if ( fParNew[i] ) delete fParNew[i];
			fParNew[i] = new TVectorD(newPar);
			TMatrixDSym newCov = GetNewCov(i);
                        if ( fCovNew[i] ) delete fCovNew[i];
			fCovNew[i] = new TMatrixDSym(newCov);
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
	//fRstart = TMath::Sqrt(fXv(0)*fXv(0) + fXv(1)*fXv(1));	// Store fit 
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
// Derivative of phases wrt initial track arameters
//
TVectorD VertexFit::DsiDa0k(Int_t i, Int_t k)
{
	// Unit 3x3 matrix
	TMatrixD M3(3, 3);
	TMatrixD Ui3(TMatrixD::kUnit, M3);
	if(i != k) Ui3.Zero();	
	//
	// Initialize D^{-1}
	TMatrixDSym D(3);	D.Zero();
	TMatrixDSym Dm1(3);
	for (Int_t k = 0; k < fNtr; k++) D += *fDi[k];
	// 
	// if vertex constraint
	if(fVtxCst) D += fCovCstInv;
	Dm1 = RegInv(D);
	// Other input variables
	TVectorD ai   = *fai[i];
	Double_t a2i  = fa2i[i];
	TMatrixDSym Wi = *fWi[i];
	TMatrixD Akt = *fAti[k];
	TMatrixD Dk  = *fDi[k];
	// final formula
	TMatrixD T = (Akt*(Dk*Dm1-Ui3))*Wi;
	TVectorD Sik = T*ai;
	Sik *= 1./a2i;
	//
	return Sik;
}
//
// Correlation matrix of new track parameters
TMatrixD VertexFit::DaiDa0k(Int_t i, Int_t k)
{
	TMatrixD M3(3, 3);
	TMatrixD M5(5, 5);
	//
	// Initialize D^{-1}
	TMatrixDSym D(3);	D.Zero();
	TMatrixDSym Dm1(3);
	for (Int_t k = 0; k < fNtr; k++) D += *fDi[k];
	// 
	// if vertex constraint
	if(fVtxCst) D += fCovCstInv;
	Dm1 = RegInv(D);
	// Other useful matrices
	TMatrixD Ait = *fAti[i];
	TMatrixD Ai(TMatrixD::kTransposed, Ait);
	//
	TMatrixD Akt = *fAti[k];
	TMatrixD Ak(TMatrixD::kTransposed, Akt);
	// i
	TMatrixD Ui3(TMatrixD::kUnit, M3);
	TMatrixD Ui5(TMatrixD::kUnit, M5);
	if (k != i) {
		Ui3.Zero();
		Ui5.Zero();
	}
	TMatrixD Mi0 = (*fDi[i]) * (Ui3 - (Dm1 * (*fDi[k])));
	TMatrixD Mik = Ait * (Mi0 * Ak);
	TMatrixD Mi = Ui5 - (*fCov[i]) * Mik;
	//
	return Mi;
}
TMatrixD VertexFit::GetNewCov(Int_t i, Int_t j)
{
	TMatrixD Cij(5,5); Cij.Zero();
	//
	// Main computation
	for(Int_t k=0; k<fNtr; k++){
		TMatrixD Mi = DaiDa0k(i, k);
		TMatrixD Mj = DaiDa0k(j, k);
		TMatrixD Mjt(TMatrixD::kTransposed,Mj);
		Cij += Mi*((*fCov[k])*Mjt);
	}
	//
	// If vertex constraint
	if(fVtxCst){
		//
		// Initialize D^{-1}
		TMatrixDSym D(3);	D.Zero();
		TMatrixDSym Dm1(3);
		for (Int_t k = 0; k < fNtr; k++) D += *fDi[k];
		D += fCovCstInv;
		Dm1 = RegInv(D);
		TMatrixD Fi = (*fCov[i])*((*fAti[i])*((*fDi[i])*Dm1));
		TMatrixD Fj = (*fCov[j])*((*fAti[j])*((*fDi[j])*Dm1));
		TMatrixD Fjt(TMatrixD::kTransposed,Fj);
		Cij += Fi*(fCovCstInv*Fjt);
	}
	//
	return Cij;
}
//
// Just diagonal terms
TMatrixDSym VertexFit::GetNewCov(Int_t i)
{
	TMatrixD  Cov = GetNewCov(i,i);
	TMatrixDSym CovSym(5);
	for(Int_t k1=0; k1<5; k1++){
		for(Int_t k2=0; k2<5; k2++)CovSym(k1,k2) = 0.5*(Cov(k1,k2)+Cov(k2,k1));
	}
	//
	return CovSym;
}
//
// Correlation parameters vertex
TMatrixD VertexFit::GetNewCovXvPar(Int_t i)
{
	TMatrixD Cxp(3,5); Cxp.Zero();
	TMatrixD M3(3,3); 
	TMatrixD M5(5,5);
	//
	// Initialize D^{-1}
	TMatrixDSym D(3);	D.Zero();
	TMatrixDSym Dm1(3);
	for(Int_t k=0; k<fNtr; k++) D += *fDi[k];
	if(fVtxCst) D += fCovCstInv;
	Dm1 = RegInv(D);
	// Other useful matrices
	//
	// Main computation
	for(Int_t k=0; k<fNtr; k++){	
/*	
		TMatrixD Akt = *fAti[k];
		TMatrixD Ak(TMatrixD::kTransposed,Akt);
		// i
		TMatrixD Ui3(TMatrixD::kUnit,M3);
		TMatrixD Ui5(TMatrixD::kUnit,M5);
		if(k != i){
			Ui3.Zero();
			Ui5.Zero();
		}	
		TMatrixD Mi0 = (*fDi[i])*(Ui3-(Dm1*(*fDi[k])));
		TMatrixD Mik = Ait*(Mi0*Ak);
		TMatrixD Mi = Ui5-(*fCov[i])*Mik;
		TMatrixD Mit(TMatrixD::kTransposed,Mi);
		Cxp += (*fDi[k])*(Ak*((*fCov[k])*Mit));
*/
	TMatrixD Mik = DaiDa0k(i, k);
	TMatrixD Mikt(TMatrixD::kTransposed, Mik);
	//std::cout<<"Mikt:"; Mikt.Print();
	TMatrixD Akt = *fAti[i];
	TMatrixD Ak(TMatrixD::kTransposed,Akt);
	//std::cout<<"Ak:"; Ak.Print();
	TMatrixDSym C0k = *fCov[k];
	//Cxp += (*fDi[k])*(Ak*(C0k*Mikt));
	TMatrixD XvAlf0k = GetDxvDpar0(k);
	Cxp += XvAlf0k*(C0k*Mikt);
	//std::cout<<"Cxp:"; Cxp.Print();
	}
	//
	if(fVtxCst){
		TMatrixD Fi = (*fCov[i])*((*fAti[i])*(*fDi[i]));
		TMatrixD Fit(TMatrixD::kTransposed, Fi);	
		Cxp += Dm1*(fCovCstInv*Fit);
	} 
	//Cxp = Dm1*Cxp;
	return Cxp;	
}
//
// Vertex derivative wrt starting paramenters
TMatrixD VertexFit::GetDxvDpar0(Int_t i)
{
	TMatrixD dXvDa0(3, 5); dXvDa0.Zero();	// Return matrix
	//
	// Initialize D^{-1}
	TMatrixDSym D(3);	D.Zero();
	TMatrixDSym Dm1(3);
	for (Int_t k = 0; k < fNtr; k++) D += *fDi[k];
	if(fVtxCst) D += fCovCstInv;
	Dm1 = RegInv(D);
	//
	// Other useful matrix
	TMatrixD Ait = *fAti[i];
	TMatrixD Ai(TMatrixD::kTransposed, Ait);
	//
	// Calculate result
	dXvDa0 = Dm1 * (*fDi[i] * Ai);
	//
	return dXvDa0;
}
//
// Handle tracks/constraints
void VertexFit::AddVtxConstraint(TVectorD xv, TMatrixDSym cov)	// Add gaussian vertex constraint
{
	//std::cout << "VertexFit::AddVtxConstraint: Not implemented yet" << std::endl;
	fVtxCst = kTRUE;				// Vertex constraint flag
	fxCst = xv;					// Constraint value
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
	Bool_t Charged = kTRUE;
	fCharged.push_back(Charged);
	//
	// Reset previous vertex temp arrays
	ResetWrkArrays();
	ffi.clear();
	fVtxDone = kFALSE;			// Reset vertex done flag
}
//
// Adding tracks one by one with charge tag
void VertexFit::AddTrk(TVectorD *par, TMatrixDSym *Cov, Bool_t Charged)			// Add track to input list
{
	fNtr++;
	fChi2List.ResizeTo(fNtr);	// Resize chi2 array
	fPar.push_back(par);			// add new track
	fCov.push_back(Cov);
	fParNew.push_back(par);			// add new track
	fCovNew.push_back(Cov);
	fCharged.push_back(Charged);
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
	fCharged.erase(fCharged.begin() + iTrk);
	//
	// Reset previous vertex temp arrays
	ResetWrkArrays();
	ffi.clear();
	fVtxDone = kFALSE;			// Reset vertex done flag
}
