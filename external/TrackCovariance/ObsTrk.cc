#include <TMath.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TDecompChol.h>
#include <TRandom.h>
#include <iostream>
#include "SolGeom.h"
#include "SolGridCov.h"
#include "ObsTrk.h"
//
// Constructors
//
// x(3) track origin, p(3) track momentum at origin, Q charge, B magnetic field in Tesla
ObsTrk::ObsTrk(TVector3 x, TVector3 p, Double_t Q, SolGridCov *GC, SolGeom *G)
{
	fB = G->B();
	SetB(fB);
	fG = G;
	fGC = GC;
	fGenX = x;
	fGenP = p;
	fGenQ = Q;
	fEflag = kFALSE;		// Electron flag 
	fEscale = 1.;			// Electron scale
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
	//
	FillGen();
	//
	fCov = CovCalc(fGenPar);
	fCovMm = CovToMm(fCov);
	fCovACTS = CovToACTS(fObsPar, fCov);
	fCovILC = CovToILC(fCov);
	//
	fObsDone = kFALSE;
}
//
// x[3] track origin, p[3] track momentum at origin, Q charge, B magnetic field in Tesla
ObsTrk::ObsTrk(Double_t *x, Double_t *p, Double_t Q, SolGridCov* GC, SolGeom *G)
{
	fB = G->B();
	SetB(fB);
	fG = G;
	fGC = GC;
	fGenX.SetXYZ(x[0],x[1],x[2]);
	fGenP.SetXYZ(p[0],p[1],p[2]);
	fGenQ = Q;
	fEflag = kFALSE;		// Electron flag 
	fEscale = 1.;			// Electron scale
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
	//
	FillGen();
	//
	fCov = CovCalc(fGenPar);
	fCovMm = CovToMm(fCov);
	fCovACTS = CovToACTS(fObsPar, fCov);
	fCovILC = CovToILC(fCov);
	//
	fObsDone = kFALSE;
}
void ObsTrk::FillGen()
{
// Fill Generated track arrays
//
	fGenPar = XPtoPar(fGenX, fGenP, fGenQ);
	fGenParMm = ParToMm(fGenPar);
	fGenParACTS = ParToACTS(fGenPar);
	fGenParILC = ParToILC(fGenPar);
}
//
void ObsTrk::FillObs()
{
// Fill Observed track arrays
//
	fObsPar = TrkUtil::CovSmear(fGenPar, fCov);
	fObsParMm = ParToMm(fObsPar);
	fObsParACTS = ParToACTS(fObsPar);
	fObsParILC = ParToILC(fObsPar);
	fObsX = ParToX(fObsPar);
	fObsP = ParToP(fObsPar);
	fObsQ = ParToQ(fObsPar);
	//
	fObsDone = kTRUE;
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
// Rescale covariance if needed
//
void ObsTrk::SetScale(Double_t scale)
{
//
// Scale covariance matrix by "scale"
//
	fEscale = scale;			// scale of resolution
	if(!fEflag){ 
		fCov *= (fEscale*fEscale);	// Rescale covariance matrix
		fEflag = kTRUE;				// Scaling flag 
		// Update covariance matrix variants
		fCovMm = CovToMm(fCov);
		fCovACTS = CovToACTS(fObsPar, fCov);
		fCovILC = CovToILC(fCov);
	}
	else std::cout<<"ObsTrk::SetScale: Already called --> no action"<<std::endl;
}

//
// Calculate covariance matrix
//
TMatrixDSym ObsTrk::CovCalc(TVectorD gPar)
{
	//
	// Check ranges
	Double_t minPt = fGC->GetMinPt();
	//if (pt < minPt) std::cout << "Warning ObsTrk::GenToObsPar: pt " << pt << " is below grid range of " << minPt << std::endl;
	Double_t maxPt = fGC->GetMaxPt();
	//if (pt > maxPt) std::cout << "Warning ObsTrk::GenToObsPar: pt " << pt << " is above grid range of " << maxPt << std::endl;
	Double_t minAn = fGC->GetMinAng();
	//if (angd < minAn) std::cout << "Warning ObsTrk::GenToObsPar: angle " << angd
	//	<< " is below grid range of " << minAn << std::endl;
	Double_t maxAn = fGC->GetMaxAng();
	//if (angd > maxAn) std::cout << "Warning ObsTrk::GenToObsPar: angle " << angd
	//	<< " is above grid range of " << maxAn << std::endl;
	//
	TMatrixDSym Cov(5);
	//
	// Check if track origin is inside beampipe and betwen the first disks
	//
	Double_t Rin = fG->GetRmin();
	Double_t ZinPos = fG->GetZminPos();
	Double_t ZinNeg = fG->GetZminNeg();
	Bool_t inside = TrkUtil::IsInside(fGenX, Rin, ZinNeg, ZinPos); // Check if in inner box
	SolTrack* trk = new SolTrack(fGenX, fGenP, fG);
	Double_t Xfirst, Yfirst, Zfirst;
	Int_t iLay = trk->FirstHit(Xfirst, Yfirst, Zfirst);
	fXfirst = TVector3(Xfirst, Yfirst, Zfirst);
  //std::cout<<"obs trk: "<<Xfirst<<","<<Yfirst<<","<<Zfirst<<std::endl;

	if (inside)
	{
		//std::cout<<"ObsTrk:: inside: x= "<<fGenX(0)<<", y= "<<fGenX(1)
                //                        <<", z= "<<fGenX(2)<<std::endl;
		// Observed track parameters
		Double_t pt = fGenP.Pt();
		Double_t angd = fGenP.Theta() * 180. / TMath::Pi();
		Cov = fGC->GetCov(pt, angd);				// Track covariance
	}
	else
	{
		//std::cout<<"ObsTrk:: outside: x= "<<fGenX(0)<<", y= "<<fGenX(1)
                //                         <<", z= "<<fGenX(2)<<std::endl;
		Bool_t Res = kTRUE; Bool_t MS = kTRUE;
		trk->CovCalc(Res, MS);					// Calculate covariance matrix
		Cov = trk->Cov();
	}					// Track covariance
	delete trk;
//
	return Cov;
}

