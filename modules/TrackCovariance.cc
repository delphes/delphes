/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2020  Universite catholique de Louvain (UCLouvain), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/** \class TrackCovariance
 *
 *  Smears track parameters according to appropriate covariance matrix.
 *
 *  \authors P. Demin - UCLouvain, Louvain-la-Neuve
 *           M. Selvaggi - CERN
 *
 */
//FIXME add reference to Bedeschi-code
//FIXME make sure about units of P, X
//FIXME fix pt > 200 GeV issue and angle > 6.41

#include "classes/DelphesClasses.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"
#include "classes/DelphesModuleFactory.h"

#include <TrackCovariance/ObsTrk.h>
#include <TrackCovariance/SolGeom.h>
#include <TrackCovariance/SolGridCov.h>

#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>

using namespace std;

class TrackCovariance: public DelphesModule
{
public:
  TrackCovariance();

  void Init() override;
  void Process() override;

private:
  Double_t fBz;
  Int_t fNMinHits;

  const std::unique_ptr<DelphesFormula> fElectronScaleFactor;
  const std::unique_ptr<DelphesFormula> fMuonScaleFactor;
  const std::unique_ptr<DelphesFormula> fChargedHadronScaleFactor;

  const std::unique_ptr<SolGeom> fGeometry;
  const std::unique_ptr<SolGridCov> fCovariance;

  AcceptanceClx *fAcx{nullptr};

  std::unique_ptr<TIterator> fItInputArray; //!

  const TObjArray *fInputArray{nullptr}; //!

  TObjArray *fOutputArray{nullptr}; //!
};

//------------------------------------------------------------------------------

TrackCovariance::TrackCovariance() :
  fElectronScaleFactor(std::make_unique<DelphesFormula>()),
  fMuonScaleFactor(std::make_unique<DelphesFormula>()),
  fChargedHadronScaleFactor(std::make_unique<DelphesFormula>()),
  fGeometry(std::make_unique<SolGeom>()),
  fCovariance(std::make_unique<SolGridCov>()) {}

//------------------------------------------------------------------------------

void TrackCovariance::Init()
{
  fBz = GetDouble("Bz", 0.0);
  fGeometry->Read(GetString("DetectorGeometry", ""));
  fGeometry->SetBz(fBz);
  fNMinHits = GetInt("NMinHits", 6);

  // scale factors to apply to resolutions
  fElectronScaleFactor->Compile(GetString("ElectronScaleFactor", "1.0"));
  fMuonScaleFactor->Compile(GetString("MuonScaleFactor", "1.0"));
  fChargedHadronScaleFactor->Compile(GetString("ChargedHadronScaleFactor", "1.0"));

  // load geometry
  fCovariance->Calc(fGeometry.get());
  fCovariance->SetMinHits(fNMinHits);
  // load geometry
  fAcx = fCovariance->AccPnt();

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "TrackMerger/tracks"));
  fItInputArray.reset(fInputArray->MakeIterator());

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
}

//------------------------------------------------------------------------------

void TrackCovariance::Process()
{
  Candidate *candidate, *mother, *particle;
  Double_t mass, p, pt, q, ct;
  Double_t dd0, ddz, dphi, dct, dp, dpt, dC;
  //
  // Get cylindrical box for fast track simulation
  //
  Double_t Rin = fGeometry->GetRmin();
  Double_t ZinPos = fGeometry->GetZminPos();
  Double_t ZinNeg = fGeometry->GetZminNeg();

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {

    // converting to meters
    particle = static_cast<Candidate *>(candidate->GetCandidates()->At(0));

    // converting to meters
    const TLorentzVector &candidatePosition = particle->Position * 1e-03;
    const TLorentzVector &candidateMomentum = particle->Momentum;

    Bool_t inside = TrkUtil::IsInside(candidatePosition.Vect(), Rin, ZinNeg, ZinPos); // Check if in inner box
    Bool_t Accept = kTRUE;
    if(inside)
      Accept = fCovariance->IsAccepted(candidateMomentum.Vect());
    else
      Accept = fCovariance->IsAccepted(candidatePosition.Vect(), candidateMomentum.Vect(), fGeometry.get());
    if(!Accept) continue;

    mass = candidateMomentum.M();

    // ********************************
    // Standard implementation with grid
    //ObsTrk track(candidatePosition.Vect(), candidateMomentum.Vect(), candidate->Charge, fCovariance.get(), fGeometry.get());
    // ********************************
    //
    // *******************************
    // Kalman implementation without any grid
    // efficiency in tracking layers use can affect acceptance
    // Comment lines below within ******** and
    // uncomment above to return to standard implementation
    //
    ObsTrk track(candidatePosition.Vect(), candidateMomentum.Vect(), candidate->Charge, fCovariance.get(), fGeometry.get());
    Int_t MinMeasure = 6; // minimum number of measurements required
    if(track.GetUmeas() < MinMeasure) continue;
    //
    // *******************************

    // apply rescaling factors to resolution
    if(TMath::Abs(candidate->PID) == 11)
    {
      track.SetScale(fElectronScaleFactor->Eval(candidateMomentum.Pt(), candidateMomentum.Eta(), candidateMomentum.Phi(), candidateMomentum.E(), candidate));
    }
    else if(TMath::Abs(candidate->PID) == 13)
    {
      track.SetScale(fMuonScaleFactor->Eval(candidateMomentum.Pt(), candidateMomentum.Eta(), candidateMomentum.Phi(), candidateMomentum.E(), candidate));
    }
    else
    {
      track.SetScale(fChargedHadronScaleFactor->Eval(candidateMomentum.Pt(), candidateMomentum.Eta(), candidateMomentum.Phi(), candidateMomentum.E(), candidate));
    }

    mother = candidate;
    candidate = static_cast<Candidate *>(candidate->Clone());

    candidate->Momentum.SetVectM(track.GetObsP(), mass);

    // converting back to mm
    candidate->InitialPosition.SetXYZT(track.GetObsX().X() * 1e03, track.GetObsX().Y() * 1e03, track.GetObsX().Z() * 1e03, candidatePosition.T() * 1e03);

    // save full covariance 5x5 matrix internally (D0, phi, Curvature, dz, ctg(theta))
    candidate->TrackCovariance = track.GetCov();

    pt = candidate->Momentum.Pt();
    p = candidate->Momentum.P();
    q = track.GetObsQ();
    ct = track.GetObsPar()[4];

    candidate->Xd = track.GetObsX().X() * 1e03;
    candidate->Yd = track.GetObsX().Y() * 1e03;
    candidate->Zd = track.GetObsX().Z() * 1e03;

    candidate->XFirstHit = track.GetFirstHit().X() * 1e03;
    candidate->YFirstHit = track.GetFirstHit().Y() * 1e03;
    candidate->ZFirstHit = track.GetFirstHit().Z() * 1e03;

    candidate->D0 = track.GetObsPar()[0] * 1e03;
    candidate->Phi = track.GetObsPar()[1];

    // inverse of curvature
    candidate->C = track.GetObsPar()[2] * 1e-03;
    candidate->DZ = track.GetObsPar()[3] * 1e03;
    candidate->CtgTheta = track.GetObsPar()[4];
    candidate->P = track.GetObsP().Mag();
    candidate->PT = pt;
    candidate->Charge = q;

    dd0 = TMath::Sqrt(track.GetCov()(0, 0)) * 1e03;
    ddz = TMath::Sqrt(track.GetCov()(3, 3)) * 1e03;
    dphi = TMath::Sqrt(track.GetCov()(1, 1));
    dct = TMath::Sqrt(track.GetCov()(4, 4));
    dpt = 2 * TMath::Sqrt(track.GetCov()(2, 2)) * pt * pt / (0.2998 * fBz);
    dp = TMath::Sqrt((1. + ct * ct) * dpt * dpt + 4 * pt * pt * ct * ct * dct * dct / (1. + ct * ct) / (1. + ct * ct));
    dC = TMath::Sqrt(track.GetCov()(2, 2)) * 1e-03;

    candidate->ErrorD0 = dd0;
    candidate->ErrorDZ = ddz;
    candidate->ErrorP = dp;
    candidate->ErrorC = dC;
    candidate->ErrorCtgTheta = dct;
    candidate->ErrorPhi = dphi;
    candidate->ErrorPT = dpt;
    //candidate->TrackResolution = dpt / pt;
    candidate->TrackResolution = dp / p;

    candidate->AddCandidate(mother);

    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("TrackCovariance", TrackCovariance);
