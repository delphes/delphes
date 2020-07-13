#include <TMath.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TDecompChol.h>
#include <TRandom.h>
#include <iostream>
#include "SolGridCov.h"
#include "ObsTrk.h"
//
// Constructors
// x(3) track origin, p(3) track momentum at origin, Q charge, B magnetic field in Tesla
ObsTrk::ObsTrk(TVector3 x, TVector3 p, Double_t Q, Double_t B, SolGridCov *GC)
{
	fGC = GC;
	fGenX = x;
	fGenP = p;
	fGenQ = Q;
	fB = B;
	fGenPar.ResizeTo(5);
	fGenParACTS.ResizeTo(6);
	fGenParILC.ResizeTo(5);
	fObsPar.ResizeTo(5);
	fObsParACTS.ResizeTo(6);
	fObsParILC.ResizeTo(5);
	fCov.ResizeTo(5, 5);
	fCovACTS.ResizeTo(6, 6);
	fCovILC.ResizeTo(5, 5);
	fGenPar = XPtoPar(x,p,Q);
	fGenParACTS = ParToACTS(fGenPar);
	fGenParILC = ParToILC(fGenPar);
	/*
	std::cout << "ObsTrk::ObsTrk: fGenPar";
	for (Int_t i = 0; i < 5; i++)std::cout << fGenPar(i) << ", ";
	std::cout << std::endl;
	*/
	fObsPar = GenToObsPar(fGenPar, fGC);
	fObsParACTS = ParToACTS(fObsPar);
	fObsParILC = ParToILC(fObsPar);
	fObsX = ParToX(fObsPar);
	fObsP = ParToP(fObsPar);
	fObsQ = ParToQ(fObsPar);
	fCovACTS = CovToACTS(fCov);
	fCovILC = CovToILC(fCov);
}
//
// Destructor
ObsTrk::~ObsTrk()
{
	fGenX.Clear();
	fGenP.Clear();
	fGenPar.Clear();
	fGenParACTS.Clear();
	fObsX.Clear();
	fObsP.Clear();
	fObsPar.Clear();
	fObsParACTS.Clear();
	fCov.Clear();
	fCovACTS.Clear();
}
TVectorD ObsTrk::XPtoPar(TVector3 x, TVector3 p, Double_t Q)
{
	//
	TVectorD Par(5);
	// Transverse parameters
	Double_t a = -Q*fB*0.2998;			// Units are Tesla, GeV and meters
	Double_t pt = p.Pt();
	Double_t C = a / (2 * pt);			// Half curvature
	//std::cout << "ObsTrk::XPtoPar: fB = " << fB << ", a = " << a << ", pt = " << pt << ", C = " << C << std::endl;
	Double_t r2 = x.Perp2();
	Double_t cross = x(0)*p(1) - x(1)*p(0);
	Double_t T = TMath::Sqrt(pt*pt - 2 * a*cross + a*a*r2);
	Double_t phi0 = TMath::ATan2((p(1) - a*x(0)) / T, (p(0) + a*x(1)) / T);	// Phi0
	Double_t D;							// Impact parameter D
	if (pt < 10.0) D = (T - pt) / a;
	else D = (-2 * cross + a*r2) / (T + pt);
	//
	Par(0) = D;		// Store D
	Par(1) = phi0;	// Store phi0
	Par(2) = C;		// Store C
	//Longitudinal parameters
	Double_t B = C*TMath::Sqrt(TMath::Max(r2 - D*D,0.0) / (1 + 2 * C*D));
	Double_t st = TMath::ASin(B) / C;
	Double_t ct = p(2) / pt;
	Double_t z0 = x(2) - ct*st;
	//
	Par(3) = z0;		// Store z0
	Par(4) = ct;		// Store cot(theta)
	//
	return Par;
}
//
TVector3 ObsTrk::ParToX(TVectorD Par)
{
	Double_t D    = Par(0);
	Double_t phi0 = Par(1);
	Double_t z0   = Par(3);
	//
	TVector3 Xval;
	Xval(0) = -D*TMath::Sin(phi0); 
	Xval(1) =  D*TMath::Cos(phi0);  
	Xval(2) =  z0;
	//
	return Xval;
}
//
TVector3 ObsTrk::ParToP(TVectorD Par)
{
	Double_t C    = Par(2);
	Double_t phi0 = Par(1);
	Double_t ct   = Par(4);
	//
	TVector3 Pval;
	Double_t pt = fB*0.2998 / TMath::Abs(2 * C);
	Pval(0) = pt*TMath::Cos(phi0);
	Pval(1) = pt*TMath::Sin(phi0);
	Pval(2) = pt*ct;
	//
	return Pval;
}
//

Double_t ObsTrk::ParToQ(TVectorD Par)
{
	return TMath::Sign(1.0, -Par(2));
}
//
TVectorD ObsTrk::GenToObsPar(TVectorD gPar, SolGridCov *GC)
{
	TVector3 p = ParToP(gPar);
	Double_t pt = p.Pt();
	Double_t tanTh = 1.0 / TMath::Abs(gPar(4));
	Double_t angd = TMath::ATan(tanTh)*180. / TMath::Pi();
	//
	// Check ranges
	Double_t minPt = GC->GetMinPt ();
	if (pt < minPt) std::cout << "Warning ObsTrk::GenToObsPar: pt " << pt << " is below grid range of " << minPt << std::endl;
	Double_t maxPt = GC->GetMaxPt();
	if (pt > maxPt) std::cout << "Warning ObsTrk::GenToObsPar: pt " << pt << " is above grid range of " << maxPt << std::endl;
	Double_t minAn = GC->GetMinAng();
	if (angd < minAn) std::cout << "Warning ObsTrk::GenToObsPar: angle " << angd 
		<< " is below grid range of " << minAn << std::endl;
	Double_t maxAn = GC->GetMaxAng();
	if (angd > maxAn) std::cout << "Warning ObsTrk::GenToObsPar: angle " << angd
		<< " is above grid range of " << maxAn << std::endl;
	//
	TMatrixDSym Cov = GC->GetCov(pt, angd);
	fCov = Cov;
	//
	// Now do Choleski decomposition and random number extraction, with appropriate stabilization
	//
	TMatrixDSym CvN = Cov;
	TMatrixDSym DCv(5); DCv.Zero();
	TMatrixDSym DCvInv(5); DCvInv.Zero();
	for (Int_t id = 0; id < 5; id++)
	{
		Double_t dVal = TMath::Sqrt(Cov(id, id));
		DCv   (id, id) = dVal;
		DCvInv(id, id) = 1.0 / dVal;
	}
	CvN.Similarity(DCvInv);			// Normalize diagonal to 1
	TDecompChol Chl(CvN);
	Bool_t OK = Chl.Decompose();		// Choleski decomposition of normalized matrix
	TMatrixD U = Chl.GetU();			// Get Upper triangular matrix
	TMatrixD Ut(TMatrixD::kTransposed, U); // Transposed of U (lower triangular)
	TVectorD r(5);
	for (Int_t i = 0; i < 5; i++)r(i) = gRandom->Gaus(0.0, 1.0);		// Array of normal random numbers
	TVectorD oPar = gPar + DCv*(Ut*r);	// Observed parameter vector
	//
	return oPar;
}
// Parameter conversion to ACTS format
TVectorD ObsTrk::ParToACTS(TVectorD Par)
{
	TVectorD pACTS(6);	// Return vector
	//
	Double_t b = -0.29988*fB / 2.;
	pACTS(0) = 1000*Par(0);		// D from m to mm
	pACTS(1) = 1000 * Par(3);	// z0 from m to mm
	pACTS(2) = Par(1);			// Phi0 is unchanged
	pACTS(3) = TMath::ATan(1.0 / Par(4)) + TMath::PiOver2();		// Theta in [0, pi] range
	pACTS(4) = Par(2) / (b*TMath::Sqrt(1 + Par(4)*Par(4)));		// q/p in GeV
	pACTS(5) = 0.0;				// Time: currently undefined
	//
	return pACTS;
}
// Covariance conversion to ACTS format
TMatrixDSym ObsTrk::CovToACTS(TMatrixDSym Cov)
{
	TMatrixDSym cACTS(6); cACTS.Zero();
	Double_t b = -0.29988*fB / 2.;
	//
	// Fill derivative matrix
	TMatrixD A(5, 5);	A.Zero();
	Double_t ct = fGenPar(4);	// cot(theta)
	Double_t C = fGenPar(2);		// half curvature
	A(0, 0) = 1000.;		// D-D	conversion to mm
	A(1, 2) = 1.0;		// phi0-phi0
	A(2, 4) = 1.0/(TMath::Sqrt(1.0 + ct*ct) * b);	// q/p-C
	A(3, 1) = 1000.;		// z0-z0 conversion to mm
	A(4, 3) = -1.0 / (1.0 + ct*ct); // theta - cot(theta)
	A(4, 4) = -C*ct / (b*pow(1.0 + ct*ct,3.0/2.0)); // q/p-cot(theta)
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

// Parameter conversion to ILC format
TVectorD ObsTrk::ParToILC(TVectorD Par)
{
	TVectorD pILC(5);	// Return vector
	//
	pILC(0) = Par(0)*1.0e3;			// d0 in mm
	pILC(1) = Par(1);				// phi0 is unchanged
	pILC(2) = -2 * Par(2)*1.0e-3;	// w in mm^-1
	pILC(3) = Par(3)*1.0e3;			// z0 in mm
	pILC(4) = Par(4);				// tan(lambda) = cot(theta)
	//
	return pILC;
}
// Covariance conversion to ILC format
TMatrixDSym ObsTrk::CovToILC(TMatrixDSym Cov)
{
	TMatrixDSym cILC(5); cILC.Zero();
	//
	// Fill derivative matrix
	TMatrixD A(5, 5);	A.Zero();
	//
	A(0, 0) = 1.0e3;		// D-d0 in mm
	A(1, 1) = 1.0;		// phi0-phi0
	A(2, 2) = -2.0e-3;	// w-C
	A(3, 3) = 1.0e3;		// z0-z0 conversion to mm
	A(4, 4) = 1.0;		// tan(lambda) - cot(theta)
	//
	TMatrixDSym Cv = Cov;
	TMatrixD At(5, 5);
	At.Transpose(A);
	Cv.Similarity(At);
	cILC = Cv;
	//
	return cILC;
}



	
