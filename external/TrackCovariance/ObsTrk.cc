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
//
// x(3) track origin, p(3) track momentum at origin, Q charge, B magnetic field in Tesla
ObsTrk::ObsTrk(TVector3 x, TVector3 p, Double_t Q, Double_t B, SolGridCov *GC)
{
	SetBfield(B);
	fGC = GC;
	fGenX = x;
	fGenP = p;
	fGenQ = Q;
	fB = B;
	fGenPar.ResizeTo(5);
	fGenParMm.ResizeTo(5);
	fGenParACTS.ResizeTo(6);
	fGenParILC.ResizeTo(5);
	fObsPar.ResizeTo(5);
	fObsParMm.ResizeTo(5);
	fObsParACTS.ResizeTo(6);
	fObsParILC.ResizeTo(5);
	fCov.ResizeTo(5, 5);
	fCovMm.ResizeTo(5, 5);
	fCovACTS.ResizeTo(6, 6);
	fCovILC.ResizeTo(5, 5);
	fGenPar = XPtoPar(x,p,Q);
	fGenParMm = ParToMm(fGenPar);
	fGenParACTS = ParToACTS(fGenPar);
	fGenParILC = ParToILC(fGenPar);
	/*
	cout << "ObsTrk::ObsTrk: fGenPar";
	for (Int_t i = 0; i < 5; i++)cout << fGenPar(i) << ", ";
	cout << endl;
	*/
	fObsPar = GenToObsPar(fGenPar, fGC);
	fObsParMm = ParToMm(fObsPar);
	fObsParACTS = ParToACTS(fObsPar);
	fObsParILC = ParToILC(fObsPar);
	fObsX = ParToX(fObsPar);
	fObsP = ParToP(fObsPar);
	fObsQ = ParToQ(fObsPar);
	fCovMm = CovToMm(fCov);
	fCovACTS = CovToACTS(fObsPar, fCov);
	fCovILC = CovToILC(fCov);
}
//
// x[3] track origin, p[3] track momentum at origin, Q charge, B magnetic field in Tesla
ObsTrk::ObsTrk(Double_t *x, Double_t *p, Double_t Q, Double_t B, SolGridCov* GC)
{
	SetBfield(B);
	fGC = GC;
	fGenX.SetXYZ(x[0],x[1],x[2]);
	fGenP.SetXYZ(p[0],p[1],p[2]);
	fGenQ = Q;
	fB = B;
	fGenPar.ResizeTo(5);
	fGenParMm.ResizeTo(5);
	fGenParACTS.ResizeTo(6);
	fGenParILC.ResizeTo(5);
	fObsPar.ResizeTo(5);
	fObsParMm.ResizeTo(5);
	fObsParACTS.ResizeTo(6);
	fObsParILC.ResizeTo(5);
	fCov.ResizeTo(5, 5);
	fCovMm.ResizeTo(5, 5);
	fCovACTS.ResizeTo(6, 6);
	fCovILC.ResizeTo(5, 5);
	fGenPar = XPtoPar(fGenX, fGenP, Q);
	fGenParACTS = ParToACTS(fGenPar);
	fGenParILC = ParToILC(fGenPar);
	/*
	cout << "ObsTrk::ObsTrk: fGenPar";
	for (Int_t i = 0; i < 5; i++)cout << fGenPar(i) << ", ";
	cout << endl;
	*/
	fObsPar = GenToObsPar(fGenPar, fGC);
	fObsParACTS = ParToACTS(fObsPar);
	fObsParILC = ParToILC(fObsPar);
	fObsX = ParToX(fObsPar);
	fObsP = ParToP(fObsPar);
	fObsQ = ParToQ(fObsPar);
	fCovACTS = CovToACTS(fObsPar, fCov);
	fCovILC = CovToILC(fCov);
}
//
// Destructor
ObsTrk::~ObsTrk()
{
	fGenX.Clear();
	fGenP.Clear();
	fGenPar.Clear();
	fGenParMm.Clear();
	fGenParACTS.Clear();
	fGenParILC.Clear();
	fObsX.Clear();
	fObsP.Clear();
	fObsPar.Clear();
	fObsParMm.Clear();
	fObsParACTS.Clear();
	fObsParILC.Clear();
	fCov.Clear();
	fCovMm.Clear();
	fCovACTS.Clear();
	fCovILC.Clear();
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
	//if (pt < minPt) cout << "Warning ObsTrk::GenToObsPar: pt " << pt << " is below grid range of " << minPt << endl;
	Double_t maxPt = GC->GetMaxPt();
	//if (pt > maxPt) cout << "Warning ObsTrk::GenToObsPar: pt " << pt << " is above grid range of " << maxPt << endl;
	Double_t minAn = GC->GetMinAng();
	//if (angd < minAn) cout << "Warning ObsTrk::GenToObsPar: angle " << angd 
	//	<< " is below grid range of " << minAn << endl;
	Double_t maxAn = GC->GetMaxAng();
	//if (angd > maxAn) cout << "Warning ObsTrk::GenToObsPar: angle " << angd
	//	<< " is above grid range of " << maxAn << endl;
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
