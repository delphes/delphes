/** \class TrackSmearing
 *
 *  Performs d0, dZ, p, Theta, Phi smearing of tracks.
 *
 *
 *
 *  \author A. Hart, M. Selvaggi
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"

#include <TFile.h>
#include <TLorentzVector.h>
#include <TProfile2D.h>
#include <TRandom3.h>

using namespace std;

class TrackSmearing: public DelphesModule
{
public:
  explicit TrackSmearing(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fBz(Steer<double>("Bz", 0.0)),
    fApplyToPileUp(Steer<bool>("ApplyToPileUp", true)),
    fD0Formula(std::make_unique<DelphesFormula>()),
    fDZFormula(std::make_unique<DelphesFormula>()),
    fPFormula(std::make_unique<DelphesFormula>()),
    fCtgThetaFormula(std::make_unique<DelphesFormula>()),
    fPhiFormula(std::make_unique<DelphesFormula>())
  {
    // read resolution formula
    // TODO: preload the files to avoid unnecessary I/O

    // !!! IF WE WANT TO KEEP ROOT INPUT !!!
    if(const std::string d0ResolutionFormula = Steer<std::string>("D0ResolutionFormula", "0.0");
      d0ResolutionFormula != "0.0")
    {
      fD0Formula->Compile(d0ResolutionFormula);
      fUseD0Formula = true;
    }
    else
    {
      fD0ResolutionFile = Steer<std::string>("D0ResolutionFile", "errors.root");
      fD0ResolutionHist = Steer<std::string>("D0ResolutionHist", "d0");
      fUseD0Formula = false;
    }

    if(const std::string dzResolutionFormula = Steer<std::string>("DZResolutionFormula", "0.0");
      dzResolutionFormula != "0.0")
    {
      fDZFormula->Compile(dzResolutionFormula);
      fUseDZFormula = true;
    }
    else
    {
      fDZResolutionFile = Steer<std::string>("DZResolutionFile", "errors.root");
      fDZResolutionHist = Steer<std::string>("DZResolutionHist", "dz");
      fUseDZFormula = false;
    }

    if(const std::string pResolutionFormula = Steer<std::string>("PResolutionFormula", "0.0");
      pResolutionFormula != "0.0")
    {
      fPFormula->Compile(pResolutionFormula);
      fUsePFormula = true;
    }
    else
    {
      fPResolutionFile = Steer<std::string>("PResolutionFile", "errors.root");
      fPResolutionHist = Steer<std::string>("PResolutionHist", "p");
      fUsePFormula = false;
    }

    if(const std::string ctgThetaFormula = Steer<std::string>("CtgThetaResolutionFormula", "0.0");
      ctgThetaFormula != "0.0")
    {
      fCtgThetaFormula->Compile(ctgThetaFormula);
      fUseCtgThetaFormula = true;
    }
    else
    {
      fCtgThetaResolutionFile = Steer<std::string>("CtgThetaResolutionFile", "errors.root");
      fCtgThetaResolutionHist = Steer<std::string>("CtgThetaResolutionHist", "ctgTheta");
      fUseCtgThetaFormula = false;
    }

    if(const std::string phiResolutionFormula = Steer<std::string>("PhiResolutionFormula", "0.0");
      phiResolutionFormula != "0.0")
    {
      fPhiFormula->Compile(phiResolutionFormula);
      fUsePhiFormula = true;
    }
    else
    {
      fPhiResolutionFile = Steer<std::string>("PhiResolutionFile", "errors.root");
      fPhiResolutionHist = Steer<std::string>("PhiResolutionHist", "phi");
      fUsePhiFormula = false;
    }
  }

  void Init() override
  {
    fInputArray = ImportArray(Steer<std::string>("InputArray", "ParticlePropagator/stableParticles"));
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "stableParticles"));
    // import beamspot
    try
    {
      fBeamSpotInputArray = ImportArray(Steer<std::string>("BeamSpotInputArray", "BeamSpotFilter/beamSpotParticle"));
    }
    catch(const std::runtime_error &)
    {
    }
  }
  void Process() override;

private:
  double ptError(const double, const double, const double, const double);

  const double fBz;
  const bool fApplyToPileUp;

  const std::unique_ptr<DelphesFormula> fD0Formula; //!
  const std::unique_ptr<DelphesFormula> fDZFormula; //!
  const std::unique_ptr<DelphesFormula> fPFormula; //!
  const std::unique_ptr<DelphesFormula> fCtgThetaFormula; //!
  const std::unique_ptr<DelphesFormula> fPhiFormula; //!

  std::string fD0ResolutionFile;
  std::string fD0ResolutionHist;
  bool fUseD0Formula;

  std::string fDZResolutionFile;
  std::string fDZResolutionHist;
  bool fUseDZFormula;

  std::string fPResolutionFile;
  std::string fPResolutionHist;
  bool fUsePFormula;

  std::string fCtgThetaResolutionFile;
  std::string fCtgThetaResolutionHist;
  bool fUseCtgThetaFormula;

  std::string fPhiResolutionFile;
  std::string fPhiResolutionHist;
  bool fUsePhiFormula;

  CandidatesCollection fInputArray; //!
  CandidatesCollection fBeamSpotInputArray; //!

  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void TrackSmearing::Process()
{
  fOutputArray->clear();

  TLorentzVector beamSpotPosition;
  double pt, eta, e, m, d0, d0Error, trueD0, dz, dzError, trueDZ, p, pError, trueP, ctgTheta, ctgThetaError, trueCtgTheta, phi, phiError, truePhi;
  double x, y, z, t, px, py, pz, theta;
  double q, r;
  double x_c, y_c, r_c, phi_0;
  double rcu, rc2, xd, yd, zd;
  const double c_light = 2.99792458E8;
  TProfile2D *d0ErrorHist = NULL,
             *dzErrorHist = NULL,
             *pErrorHist = NULL,
             *ctgThetaErrorHist = NULL,
             *phiErrorHist = NULL;

  if(fBeamSpotInputArray->empty())
    beamSpotPosition.SetXYZT(0.0, 0.0, 0.0, 0.0);
  else
  {
    Candidate &beamSpotCandidate = *((Candidate *)fBeamSpotInputArray->at(0));
    beamSpotPosition = beamSpotCandidate.Position;
  }

  if(!fUseD0Formula)
  {
    TFile *fin = TFile::Open(fD0ResolutionFile.c_str());
    d0ErrorHist = (TProfile2D *)fin->Get(fD0ResolutionHist.c_str());
    d0ErrorHist->SetDirectory(0);
    fin->Close();
  }
  if(!fUseDZFormula)
  {
    TFile *fin = TFile::Open(fDZResolutionFile.c_str());
    dzErrorHist = (TProfile2D *)fin->Get(fDZResolutionHist.c_str());
    dzErrorHist->SetDirectory(0);
    fin->Close();
  }
  if(!fUsePFormula)
  {
    TFile *fin = TFile::Open(fPResolutionFile.c_str());
    pErrorHist = (TProfile2D *)fin->Get(fPResolutionHist.c_str());
    pErrorHist->SetDirectory(0);
    fin->Close();
  }
  if(!fUseCtgThetaFormula)
  {
    TFile *fin = TFile::Open(fCtgThetaResolutionFile.c_str());
    ctgThetaErrorHist = (TProfile2D *)fin->Get(fCtgThetaResolutionHist.c_str());
    ctgThetaErrorHist->SetDirectory(0);
    fin->Close();
  }
  if(!fUsePhiFormula)
  {
    TFile *fin = TFile::Open(fPhiResolutionFile.c_str());
    phiErrorHist = (TProfile2D *)fin->Get(fPhiResolutionHist.c_str());
    phiErrorHist->SetDirectory(0);
    fin->Close();
  }

  for(Candidate *const &candidate : *fInputArray)
  {
    const TLorentzVector &momentum = candidate->Momentum;
    const TLorentzVector &position = candidate->InitialPosition;

    pt = momentum.Pt();
    eta = momentum.Eta();
    e = momentum.E();
    m = momentum.M();

    d0 = trueD0 = candidate->D0;
    dz = trueDZ = candidate->DZ;

    p = trueP = candidate->P;
    ctgTheta = trueCtgTheta = candidate->CtgTheta;
    phi = truePhi = candidate->Phi;

    if(fUseD0Formula)
      d0Error = fD0Formula->Eval(pt, eta, phi, e, candidate);
    else
    {
      int xbin, ybin;

      xbin = pt < d0ErrorHist->GetXaxis()->GetXmax() ? d0ErrorHist->GetXaxis()->FindBin(pt) : d0ErrorHist->GetXaxis()->GetBinCenter(d0ErrorHist->GetXaxis()->GetNbins());
      ybin = d0ErrorHist->GetYaxis()->FindBin(std::fabs(eta));
      d0Error = d0ErrorHist->GetBinContent(xbin, ybin);
      if(!d0Error)
        d0Error = -1.0;
    }
    if(d0Error < 0.0)
      continue;

    if(fUseDZFormula)
      dzError = fDZFormula->Eval(pt, eta, phi, e, candidate);
    else
    {
      int xbin, ybin;

      xbin = pt < dzErrorHist->GetXaxis()->GetXmax() ? dzErrorHist->GetXaxis()->FindBin(pt) : dzErrorHist->GetXaxis()->GetBinCenter(dzErrorHist->GetXaxis()->GetNbins());
      ybin = dzErrorHist->GetYaxis()->FindBin(std::fabs(eta));
      dzError = dzErrorHist->GetBinContent(xbin, ybin);
      if(!dzError)
        dzError = -1.0;
    }
    if(dzError < 0.0)
      continue;

    if(fUsePFormula)
      pError = fPFormula->Eval(pt, eta, phi, e, candidate) * p;
    else
    {
      int xbin, ybin;

      xbin = pt < pErrorHist->GetXaxis()->GetXmax() ? pErrorHist->GetXaxis()->FindBin(pt) : pErrorHist->GetXaxis()->GetBinCenter(pErrorHist->GetXaxis()->GetNbins());
      ybin = pErrorHist->GetYaxis()->FindBin(std::fabs(eta));
      pError = pErrorHist->GetBinContent(xbin, ybin) * p;
      if(!pError)
        pError = -1.0;
    }
    if(pError < 0.0)
      continue;

    if(fUseCtgThetaFormula)
      ctgThetaError = fCtgThetaFormula->Eval(pt, eta, phi, e, candidate);
    else
    {
      int xbin, ybin;

      xbin = pt < ctgThetaErrorHist->GetXaxis()->GetXmax() ? ctgThetaErrorHist->GetXaxis()->FindBin(pt) : ctgThetaErrorHist->GetXaxis()->GetBinCenter(ctgThetaErrorHist->GetXaxis()->GetNbins());
      ybin = ctgThetaErrorHist->GetYaxis()->FindBin(std::fabs(eta));
      ctgThetaError = ctgThetaErrorHist->GetBinContent(xbin, ybin);
      if(!ctgThetaError)
        ctgThetaError = -1.0;
    }
    if(ctgThetaError < 0.0)
      continue;

    if(fUsePhiFormula)
      phiError = fPhiFormula->Eval(pt, eta, phi, e, candidate);
    else
    {
      int xbin, ybin;

      xbin = pt < phiErrorHist->GetXaxis()->GetXmax() ? phiErrorHist->GetXaxis()->FindBin(pt) : phiErrorHist->GetXaxis()->GetBinCenter(phiErrorHist->GetXaxis()->GetNbins());
      ybin = phiErrorHist->GetYaxis()->FindBin(std::fabs(eta));
      phiError = phiErrorHist->GetBinContent(xbin, ybin);
      if(!phiError)
        phiError = -1.0;
    }
    if(phiError < 0.0)
      continue;

    if(fApplyToPileUp || !candidate->IsPU)
    {
      d0 = gRandom->Gaus(d0, d0Error);
      dz = gRandom->Gaus(dz, dzError);
      p = gRandom->Gaus(p, pError);
      ctgTheta = gRandom->Gaus(ctgTheta, ctgThetaError);
      phi = gRandom->Gaus(phi, phiError);
    }

    if(p < 0.0) continue;
    while(phi > M_PI) phi -= 2. * M_PI;
    while(phi <= -M_PI) phi += 2. * M_PI;

    Candidate *new_candidate = static_cast<Candidate *>(candidate->Clone());
    new_candidate->D0 = d0;
    new_candidate->DZ = dz;
    new_candidate->P = p;
    new_candidate->CtgTheta = ctgTheta;
    new_candidate->Phi = phi;

    theta = std::acos(ctgTheta / std::sqrt(1.0 + ctgTheta * ctgTheta));
    new_candidate->Momentum.SetPx(p * std::cos(phi) * std::sin(theta));
    new_candidate->Momentum.SetPy(p * std::sin(phi) * std::sin(theta));
    new_candidate->Momentum.SetPz(p * std::cos(theta));
    new_candidate->Momentum.SetE(std::sqrt(p * p + m * m));
    new_candidate->PT = new_candidate->Momentum.Pt();

    x = position.X();
    y = position.Y();
    z = position.Z();
    t = position.T();
    px = new_candidate->Momentum.Px();
    py = new_candidate->Momentum.Py();
    pz = new_candidate->Momentum.Pz();
    pt = new_candidate->Momentum.Pt();

    // -- solve for delta: d0' = ( (x+delta)*py' - (y+delta)*px' )/pt'

    new_candidate->InitialPosition.SetX(x + ((px * y - py * x + d0 * pt) / (py - px)));
    new_candidate->InitialPosition.SetY(y + ((px * y - py * x + d0 * pt) / (py - px)));
    x = new_candidate->InitialPosition.X();
    y = new_candidate->InitialPosition.Y();
    new_candidate->InitialPosition.SetZ(z + ((pz * (px * (x - beamSpotPosition.X()) + py * (y - beamSpotPosition.Y())) + pt * pt * (dz - z)) / (pt * pt)));
    z = new_candidate->InitialPosition.Z();

    new_candidate->InitialPosition.SetT(t);

    // update closest approach
    x *= 1.0E-3;
    y *= 1.0E-3;
    z *= 1.0E-3;

    q = new_candidate->Charge;

    r = pt / (q * fBz) * 1.0E9 / c_light; // in [m]
    phi_0 = std::atan2(py, px); // [rad] in [-pi, pi]

    // 2. helix axis coordinates
    x_c = x + r * std::sin(phi_0);
    y_c = y - r * std::cos(phi_0);
    r_c = std::hypot(x_c, y_c);

    rcu = std::fabs(r);
    rc2 = r_c * r_c;

    // calculate coordinates of closest approach to track circle in transverse plane xd, yd, zd
    xd = x_c * x_c * x_c - x_c * rcu * r_c + x_c * y_c * y_c;
    xd = (rc2 > 0.0) ? xd / rc2 : -999;
    yd = y_c * (-rcu * r_c + rc2);
    yd = (rc2 > 0.0) ? yd / rc2 : -999;
    zd = z + (std::sqrt(xd * xd + yd * yd) - std::sqrt(x * x + y * y)) * pz / pt;

    new_candidate->Xd = xd * 1.0E3;
    new_candidate->Yd = yd * 1.0E3;
    new_candidate->Zd = zd * 1.0E3;

    if(fApplyToPileUp || !new_candidate->IsPU)
    {
      new_candidate->ErrorD0 = d0Error;
      new_candidate->ErrorDZ = dzError;
      new_candidate->ErrorP = pError;
      new_candidate->ErrorCtgTheta = ctgThetaError;
      new_candidate->ErrorPhi = phiError;
      new_candidate->ErrorPT = ptError(p, ctgTheta, pError, ctgThetaError);
      new_candidate->TrackResolution = pError / p;
    }

    new_candidate->AddCandidate(candidate);
    fOutputArray->emplace_back(new_candidate);
  }
}

double TrackSmearing::ptError(const double p, const double ctgTheta, const double dP, const double dCtgTheta)
{
  double a, b;
  a = (p * p * ctgTheta * ctgTheta * dCtgTheta * dCtgTheta) / ((ctgTheta * ctgTheta + 1) * (ctgTheta * ctgTheta + 1) * (ctgTheta * ctgTheta + 1));
  b = (dP * dP) / (ctgTheta * ctgTheta + 1);
  return sqrt(a + b);
}

//------------------------------------------------------------------------------

REGISTER_MODULE("TrackSmearing", TrackSmearing);
