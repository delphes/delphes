#include <iostream>

#include <TMath.h>
#include <TVectorD.h>
#include <TVector3.h>
#include <TMatrixDSym.h>
#include <TDecompChol.h>
#include <TMatrixDSymEigen.h>

#include "SolGridCov.h"
#include "SolGeom.h"
#include "SolTrack.h"

using namespace std;

SolGridCov::SolGridCov()
{
  // Define pt-polar angle grid
  fNpt = 22;
  fPta.ResizeTo(fNpt);
  Double_t p[] = { 0.1, 0.2, 0.5, 0.7, 1., 2., 3., 4., 6., 8., 10., 15.,
                   20., 25., 30., 40., 50., 60., 80., 100., 150., 200. };
  for (Int_t ip = 0; ip < fNpt; ip++) fPta(ip) = p[ip];

  fNang = 13;
  fAnga.ResizeTo(fNang);
  Double_t a[] = { 10., 15., 20., 25., 30., 35., 40., 45., 50., 60., 70., 80., 90. };
  for (Int_t ia = 0; ia < fNang; ia++) fAnga(ia) = a[ia];
  fCov = new TMatrixDSym[fNpt * fNang];
  for (Int_t ip = 0; ip < fNpt; ip++)
  {
    for (Int_t ia = 0; ia < fNang; ia++) fCov[ip * fNang + ia].ResizeTo(5, 5);
  }
}

SolGridCov::~SolGridCov()
{
  delete[] fCov;
  delete fAcc;
}

void SolGridCov::Calc(SolGeom *G)
{
  TVectorD pta = fPta;
  TVectorD anga = fAnga;
  Bool_t Res = kTRUE; Bool_t MS = kTRUE; // Resolution and multiple scattering flags
  for (Int_t ip = 0; ip < fNpt; ip++) // Loop on pt grid
  {
    Int_t ipt = TMath::Nint(10 * pta(ip));
    for (Int_t ia = 0; ia < fNang; ia++) // Loop on angle grid
    {
      Double_t th = TMath::Pi() * (anga(ia)) / 180.;
      Double_t x[3], p[3];
      x[0] = 0; x[1] = 0; x[2] = 0; // Set origin
      p[0] = pta(ip); p[1] = 0; p[2] = pta(ip) / TMath::Tan(th);
      //
      SolTrack *tr = new SolTrack(x, p, G); // Initialize track
      tr->CovCalc(Res, MS); // Calculate covariance
      fCov[ip * fNang + ia] = tr->Cov(); // Get covariance
    }
  }

// Now make acceptance
fAcc = new AcceptanceClx(G);
}


//
Bool_t SolGridCov::IsAccepted(Double_t pt, Double_t Theta)
{
	//
	// pt in GeV, Theta in degrees
	//
	Bool_t Accept = kFALSE;
	if (fAcc->HitNumber(pt, Theta) >= fNminHits)Accept = kTRUE;
	//
	return Accept;
}
//
Bool_t SolGridCov::IsAccepted(Double_t *p)
{
	//
	// pt in GeV, Theta in degrees
	//
	Bool_t Accept = kFALSE;
	Double_t pt = TMath::Sqrt(p[0] * p[0] + p[1] * p[1]);
	Double_t th = 180. * TMath::ATan2(pt, p[2])/TMath::Pi();
	if (fAcc->HitNumber(pt,th) >= fNminHits)Accept = kTRUE;
	//
	return Accept;
}
//
Bool_t SolGridCov::IsAccepted(TVector3 p)
{
	//
	// pt in GeV, Theta in degrees
	//
	Bool_t Accept = kFALSE;
	Double_t pt = p.Pt();
	Double_t th = 180.*TMath::ACos(p.CosTheta())/TMath::Pi();
	if (fAcc->HitNumber(pt,th) >= fNminHits)Accept = kTRUE;
	//
	return Accept;
}
//
// Detailed acceptance check
//
Bool_t SolGridCov::IsAccepted(TVector3 x, TVector3 p, SolGeom* G)
{
	Bool_t Accept = kFALSE;
	//
	// Check if track origin is inside beampipe and betwen the first disks
	//
	Double_t Rin = G->GetRmin();
	Double_t ZinPos = G->GetZminPos();
	Double_t ZinNeg = G->GetZminNeg();
	Bool_t inside = TrkUtil::IsInside(x, Rin, ZinNeg, ZinPos); // Check if in inner box
	if (inside) Accept = IsAccepted(p);
	else
	{
		SolTrack* trk = new SolTrack(x, p, G);
		if (trk->nmHit() >= fNminHits)Accept = kTRUE;
		delete trk;
	}
	//
	return Accept;
}

//
// Find bin in grid
Int_t SolGridCov::GetMinIndex(Double_t xval, Int_t N, TVectorD x)
{
  Int_t min = -1; // default for xval below the lower limit
  if (xval < x(0))return min;
  if (xval>x(N - 1)){ min = N; return min; }
  for (Int_t i = 0; i < N; i++) if (xval>x(i))min = i;
  return min;
}
// Force positive definitness in normalized matrix
TMatrixDSym SolGridCov::MakePosDef(TMatrixDSym NormMat)
{
  // Input: symmetric matrix with 1's on diagonal
  // Output: positive definite matrix with 1's on diagonal

  // Default return value
  TMatrixDSym rMatN = NormMat;
  // Check the diagonal
  Bool_t Check = kFALSE;
  Int_t Size = NormMat.GetNcols();
  for (Int_t i = 0; i < Size; i++)if (TMath::Abs(NormMat(i, i) - 1.0)>1.0E-15)Check = kTRUE;
  if (Check)
  {
    cout << "SolGridCov::MakePosDef: input matrix doesn ot have 1 on diagonal. Abort." << endl;
    return rMatN;
  }
  // Diagonalize matrix
  TMatrixDSymEigen Eign(NormMat);
  TMatrixD U = Eign.GetEigenVectors();
  TVectorD lambda = Eign.GetEigenValues();
  // Reset negative eigenvalues to small positive value
  TMatrixDSym D(Size); D.Zero(); Double_t eps = 1.0e-13;
  for (Int_t i = 0; i < Size; i++)
  {
    D(i, i) = lambda(i);
    if (lambda(i) <= 0) D(i, i) = eps;
  }
  // Rebuild matrix
  TMatrixD Ut(TMatrixD::kTransposed, U);
  TMatrixD rMat = (U*D)*Ut;       // Now it is positive defite
  // Restore all ones on diagonal
  for (Int_t i1 = 0; i1 < Size; i1++)
  {
    Double_t rn1 = TMath::Sqrt(rMat(i1, i1));
    for (Int_t i2 = 0; i2 <= i1; i2++)
    {
      Double_t rn2 = TMath::Sqrt(rMat(i2, i2));
      rMatN(i1, i2) = 0.5*(rMat(i1, i2) + rMat(i2, i1)) / (rn1*rn2);
      rMatN(i2, i1) = rMatN(i1, i2);
    }
  }
  return rMatN;
}
// Interpolate covariance matrix: Bi-linear interpolation
TMatrixDSym SolGridCov::GetCov(Double_t pt, Double_t ang)
{
  // pt in GeV and ang in degrees
  Int_t minPt = GetMinIndex(pt, fNpt, fPta);
  if (minPt == -1)minPt = 0;
  if (minPt == fNpt - 1)minPt = fNpt - 2;
  Double_t dpt = fPta(minPt + 1) - fPta(minPt);
  // Put ang in 0-90 range
  ang = TMath::Abs(ang);
  while (ang > 90.)ang -= 90.;  // Needs to be fixed
  Int_t minAng = GetMinIndex(ang, fNang, fAnga);
  if (minAng == -1)minAng = 0;
  if (minAng == fNang - 1)minAng = fNang - 2;
  Double_t dang = fAnga(minAng + 1) - fAnga(minAng);
  //
  Double_t tpt = (pt - fPta(minPt)) / dpt;
  Double_t tang = (ang - fAnga(minAng)) / dang;
  //
  TMatrixDSym C11 = fCov[minPt * fNang + minAng];
  TMatrixDSym C12 = fCov[minPt * fNang + minAng + 1];
  TMatrixDSym C21 = fCov[(minPt + 1) * fNang + minAng];
  TMatrixDSym C22 = fCov[(minPt + 1) * fNang + minAng + 1];
  TMatrixDSym Cv = ((1-tpt) * (1-tang)) * C11 +
                   ((1-tpt) *    tang ) * C12 +
                   (   tpt  * (1-tang)) * C21 +
                   (   tpt  *    tang ) * C22;
  // Check for positive definiteness
  TMatrixDSym CvN = Cv;
  TMatrixDSym DCvInv(5); DCvInv.Zero();
  for (Int_t id = 0; id < 5; id++) DCvInv(id, id) = 1.0 / TMath::Sqrt(Cv(id, id));
  CvN.Similarity(DCvInv); // Normalize diagonal to 1
  TDecompChol Chl(CvN);
  if (!Chl.Decompose())
  {
    std::cout << "SolGridCov::GetCov: Interpolated matrix not positive definite. Recovering ...." << std::endl;
    TMatrixDSym rCv = MakePosDef(CvN); CvN = rCv;
    TMatrixDSym DCv(5); DCv.Zero();
    for (Int_t id = 0; id < 5; id++) DCv(id, id) = TMath::Sqrt(Cv(id, id));
    Cv = CvN.Similarity(DCv); // Normalize diagonal to 1
  }

  return Cv;
}
