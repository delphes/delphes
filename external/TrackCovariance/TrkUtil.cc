#include "TrkUtil.h"
#include <TMath.h>
#include <iostream>

// Constructor
TrkUtil::TrkUtil(Double_t Bz)
{
	fBz = Bz;
}
TrkUtil::TrkUtil()
{
	fBz = 0.0;
}
//
// Destructor
TrkUtil::~TrkUtil()
{
	fBz = 0.0;
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
	//cout << "ObsTrk::XPtoPar: fB = " << fB << ", a = " << a << ", pt = " << pt << ", C = " << C << endl;
	Double_t r2 = x.Perp2();
	Double_t cross = x(0) * p(1) - x(1) * p(0);
	Double_t T = TMath::Sqrt(pt * pt - 2 * a * cross + a * a * r2);
	Double_t phi0 = TMath::ATan2((p(1) - a * x(0)) / T, (p(0) + a * x(1)) / T);	// Phi0
	Double_t D;							// Impact parameter D
	if (pt < 10.0) D = (T - pt) / a;
	else D = (-2 * cross + a * r2) / (T + pt);
	//
	Par(0) = D;		// Store D
	Par(1) = phi0;	// Store phi0
	Par(2) = C;		// Store C
	//Longitudinal parameters
	Double_t B = C * TMath::Sqrt(TMath::Max(r2 - D * D, 0.0) / (1 + 2 * C * D));
	Double_t st = TMath::ASin(B) / C;
	Double_t ct = p(2) / pt;
	Double_t z0 = x(2) - ct * st;
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
	Xval(0) = -D * TMath::Sin(phi0);
	Xval(1) = D * TMath::Cos(phi0);
	Xval(2) = z0;
	//
	return Xval;
}
//
TVector3 TrkUtil::ParToP(TVectorD Par)
{
	if (fBz == 0.0)
std::cout << "TrkUtil::ParToP: Warning Bz not set" << std::endl;
	//
	return ParToP(Par,fBz);
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
	Pval(0) = pt * TMath::Cos(phi0);
	Pval(1) = pt * TMath::Sin(phi0);
	Pval(2) = pt * ct;
	//
	return Pval;
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
	pACTS(1) = 1000 * Par(3);		// z0 from m to mm
	pACTS(2) = Par(1);			// Phi0 is unchanged
	pACTS(3) = TMath::ATan2(1.0, Par(4));		// Theta in [0, pi] range
	pACTS(4) = Par(2) / (b * TMath::Sqrt(1 + Par(4) * Par(4)));		// q/p in GeV
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
	A(2, 4) = 1.0 / (TMath::Sqrt(1.0 + ct * ct) * b);	// q/p-C
	A(3, 1) = 1000.;		// z0-z0 conversion to mm
	A(4, 3) = -1.0 / (1.0 + ct * ct); // theta - cot(theta)
	A(4, 4) = -C * ct / (b * TMath::Power(1.0 + ct * ct, 3.0 / 2.0)); // q/p-cot(theta)
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
}
