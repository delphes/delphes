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

#include <TrackCovariance/ObsTrk.h>
#include <TrackCovariance/SolGeom.h>
#include <TrackCovariance/SolGridCov.h>

#include <TLorentzVector.h>

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

  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!
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

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
}

//------------------------------------------------------------------------------

void TrackCovariance::Process()
{
  fOutputArray->clear();

  Double_t mass, p, pt, q, ct;
  Double_t dd0, ddz, dphi, dct, dp, dpt, dC;
  //
  // Get cylindrical box for fast track simulation
  //
  Double_t Rin = fGeometry->GetRmin();
  Double_t ZinPos = fGeometry->GetZminPos();
  Double_t ZinNeg = fGeometry->GetZminNeg();

  for(Candidate *const &candidate : *fInputArray)
  {

    // converting to meters
    Candidate *particle = static_cast<Candidate *>(candidate->GetCandidates().at(0));

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
    if(std::abs(candidate->PID) == 11)
    {
      track.SetScale(fElectronScaleFactor->Eval(candidateMomentum.Pt(), candidateMomentum.Eta(), candidateMomentum.Phi(), candidateMomentum.E(), candidate));
    }
    else if(std::abs(candidate->PID) == 13)
    {
      track.SetScale(fMuonScaleFactor->Eval(candidateMomentum.Pt(), candidateMomentum.Eta(), candidateMomentum.Phi(), candidateMomentum.E(), candidate));
    }
    else
    {
      track.SetScale(fChargedHadronScaleFactor->Eval(candidateMomentum.Pt(), candidateMomentum.Eta(), candidateMomentum.Phi(), candidateMomentum.E(), candidate));
    }

    Candidate *new_candidate = static_cast<Candidate *>(candidate->Clone());

    new_candidate->Momentum.SetVectM(track.GetObsP(), mass);

    // converting back to mm
    new_candidate->InitialPosition.SetXYZT(track.GetObsX().X() * 1e03, track.GetObsX().Y() * 1e03, track.GetObsX().Z() * 1e03, candidatePosition.T() * 1e03);

    // save full covariance 5x5 matrix internally (D0, phi, Curvature, dz, ctg(theta))
    new_candidate->TrackCovariance = track.GetCov();

    pt = new_candidate->Momentum.Pt();
    p = new_candidate->Momentum.P();
    q = track.GetObsQ();
    ct = track.GetObsPar()[4];

    new_candidate->Xd = track.GetObsX().X() * 1e03;
    new_candidate->Yd = track.GetObsX().Y() * 1e03;
    new_candidate->Zd = track.GetObsX().Z() * 1e03;

    new_candidate->XFirstHit = track.GetFirstHit().X() * 1e03;
    new_candidate->YFirstHit = track.GetFirstHit().Y() * 1e03;
    new_candidate->ZFirstHit = track.GetFirstHit().Z() * 1e03;

    new_candidate->D0 = track.GetObsPar()[0] * 1e03;
    new_candidate->Phi = track.GetObsPar()[1];

    // inverse of curvature
    new_candidate->C = track.GetObsPar()[2] * 1e-03;
    new_candidate->DZ = track.GetObsPar()[3] * 1e03;
    new_candidate->CtgTheta = track.GetObsPar()[4];
    new_candidate->P = track.GetObsP().Mag();
    new_candidate->PT = pt;
    new_candidate->Charge = q;

    dd0 = std::sqrt(track.GetCov()(0, 0)) * 1e03;
    ddz = std::sqrt(track.GetCov()(3, 3)) * 1e03;
    dphi = std::sqrt(track.GetCov()(1, 1));
    dct = std::sqrt(track.GetCov()(4, 4));
    dpt = 2 * std::sqrt(track.GetCov()(2, 2)) * pt * pt / (0.2998 * fBz);
    dp = std::sqrt((1. + ct * ct) * dpt * dpt + 4 * pt * pt * ct * ct * dct * dct / (1. + ct * ct) / (1. + ct * ct));
    dC = std::sqrt(track.GetCov()(2, 2)) * 1e-03;

    new_candidate->ErrorD0 = dd0;
    new_candidate->ErrorDZ = ddz;
    new_candidate->ErrorP = dp;
    new_candidate->ErrorC = dC;
    new_candidate->ErrorCtgTheta = dct;
    new_candidate->ErrorPhi = dphi;
    new_candidate->ErrorPT = dpt;
    //new_candidate->TrackResolution = dpt / pt;
    new_candidate->TrackResolution = dp / p;

    new_candidate->AddCandidate(candidate);

    fOutputArray->emplace_back(new_candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("TrackCovariance", TrackCovariance);
