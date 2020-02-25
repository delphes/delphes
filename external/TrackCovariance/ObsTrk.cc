#include <iostream>

#include <TMath.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TDecompChol.h>
#include <TRandom.h>

#include "SolGridCov.h"
#include "ObsTrk.h"

using namespace std;

// x(3) track origin, p(3) track momentum at origin, Q charge, B magnetic field in Tesla
ObsTrk::ObsTrk(TVector3 x, TVector3 p, Double_t Q, Double_t B, SolGridCov *GC)
{
  fGC = GC;
  fGenX = x;
  fGenP = p;
  fGenQ = Q;
  fB = B;
  fGenPar.ResizeTo(5);
  fObsPar.ResizeTo(5);
  fCov.ResizeTo(5, 5);
  fGenPar = XPtoPar(x, p, Q);
  fObsPar = GenToObsPar(fGenPar, fGC);
  fObsX = ParToX(fObsPar);
  fObsP = ParToP(fObsPar);
  fObsQ = ParToQ(fObsPar);
}

ObsTrk::~ObsTrk()
{
}

TVectorD ObsTrk::XPtoPar(TVector3 x, TVector3 p, Double_t Q)
{
  TVectorD Par(5);
  // Transverse parameters
  Double_t a = -Q * fB * 0.2998; // Units are Tesla, GeV and m
  Double_t pt = p.Pt();
  Double_t C = a / (2 * pt); // Half curvature

  Double_t r2 = x.Perp2();
  Double_t cross = x(0) * p(1) - x(1) * p(0);
  Double_t T = TMath::Sqrt(pt * pt - 2 * a * cross + a * a * r2);
  Double_t phi0 = TMath::ATan2((p(1) - a * x(0)) / T, (p(0) + a * x(1)) / T); // Phi0
  Double_t D; // Impact parameter D
  if (pt < 10.0) D = (T - pt) / a;
  else D = (-2 * cross + a * r2) / (T + pt);

  Par(0) = D; // Store D
  Par(1) = phi0; // Store phi0
  Par(2) = C; // Store C
  // Longitudinal parameters
  Double_t B = C * TMath::Sqrt(TMath::Max(r2 - D * D,0.0) / (1 + 2 * C * D));
  Double_t st = TMath::ASin(B) / C;
  Double_t ct = p(2) / pt;
  Double_t z0 = x(2) - ct * st;

  Par(3) = z0; // Store z0
  Par(4) = ct; // Store cot(theta)

  return Par;
}

TVector3 ObsTrk::ParToX(TVectorD Par)
{
  Double_t D    = Par(0);
  Double_t phi0 = Par(1);
  Double_t z0   = Par(3);

  TVector3 Xval;
  Xval(0) = -D * TMath::Sin(phi0);
  Xval(1) =  D * TMath::Cos(phi0);
  Xval(2) =  z0;

  return Xval;
}

TVector3 ObsTrk::ParToP(TVectorD Par)
{
  Double_t C    = Par(2);
  Double_t phi0 = Par(1);
  Double_t ct   = Par(4);
  //
  TVector3 Pval;
  Double_t pt = fB * 0.2998 / TMath::Abs(2 * C);
  Pval(0) = pt * TMath::Cos(phi0);
  Pval(1) = pt * TMath::Sin(phi0);
  Pval(2) = pt * ct;
  //
  return Pval;
}

Double_t ObsTrk::ParToQ(TVectorD Par)
{
  return TMath::Sign(1.0, -Par(2));
}

TVectorD ObsTrk::GenToObsPar(TVectorD gPar, SolGridCov *GC)
{
  TVector3 p = ParToP(gPar);
  Double_t pt = p.Pt();
  Double_t tanTh = 1.0 / TMath::Abs(gPar(4));
  Double_t angd = TMath::ATan(tanTh) * 180. / TMath::Pi();
  // Check ranges
  Double_t minPt = GC->GetMinPt ();
  if (pt < minPt) cout << "Warning ObsTrk::GenToObsPar: pt " << pt << " is below grid range of " << minPt << endl;
  Double_t maxPt = GC->GetMaxPt();
  if (pt > maxPt) cout << "Warning ObsTrk::GenToObsPar: pt " << pt << " is above grid range of " << maxPt << endl;
  Double_t minAn = GC->GetMinAng();
  if (angd < minAn) cout << "Warning ObsTrk::GenToObsPar: angle " << angd
    << " is below grid range of " << minAn << endl;
  Double_t maxAn = GC->GetMaxAng();
  if (angd > maxAn) cout << "Warning ObsTrk::GenToObsPar: angle " << angd
    << " is above grid range of " << maxAn << endl;
  TMatrixDSym Cov = GC->GetCov(pt, angd);
  fCov = Cov;
  // Now do Choleski decomposition and random number extraction, with appropriate stabilization
  TMatrixDSym CvN = Cov;
  TMatrixDSym DCv(5); DCv.Zero();
  TMatrixDSym DCvInv(5); DCvInv.Zero();
  for (Int_t id = 0; id < 5; id++)
  {
    Double_t dVal = TMath::Sqrt(Cov(id, id));
    DCv   (id, id) = dVal;
    DCvInv(id, id) = 1.0 / dVal;
  }
  CvN.Similarity(DCvInv); // Normalize diagonal to 1
  TDecompChol Chl(CvN);
  Bool_t OK = Chl.Decompose(); // Choleski decomposition of normalized matrix
  TMatrixD U = Chl.GetU(); // Get Upper triangular matrix
  TMatrixD Ut(TMatrixD::kTransposed, U); // Transposed of U (lower triangular)
  TVectorD r(5);
  for (Int_t i = 0; i < 5; i++) r(i) = gRandom->Gaus(0.0, 1.0); // Array of normal random numbers
  TVectorD oPar = gPar + DCv * (Ut * r); // Observed parameter vector

  return oPar;
}
