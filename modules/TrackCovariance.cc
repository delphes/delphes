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
 *  \authors F. Bedeschi - INFN Pisa
*            P. Demin - UCLouvain, Louvain-la-Neuve
 *           M. Selvaggi - CERN
 *
 *
 */

//FIXME add reference to Bedeschi-code
//FIXME make sure about units of P, X
//FIXME fix pt > 200 GeV issue and angle > 6.41

#include "modules/TrackCovariance.h"

#include "classes/DelphesClasses.h"

#include "TrackCovariance/ObsTrk.h"
#include "TrackCovariance/SolGeom.h"
#include "TrackCovariance/SolGridCov.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "TMath.h"
#include "TObjArray.h"

#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

TrackCovariance::TrackCovariance() :
  fElectronScaleFactor(0), fMuonScaleFactor(0), fChargedHadronScaleFactor(0),
  fGeometry(0), fCovariance(0), fAcx(0)
{
  fGeometry = new SolGeom();
  fCovariance = new SolGridCov();
  fElectronScaleFactor = new DelphesFormula;
  fMuonScaleFactor = new DelphesFormula;
  fChargedHadronScaleFactor = new DelphesFormula;
}

//------------------------------------------------------------------------------

TrackCovariance::~TrackCovariance()
{
  if(fGeometry) delete fGeometry;
  if(fCovariance) delete fCovariance;
  if(fElectronScaleFactor) delete fElectronScaleFactor;
  if(fMuonScaleFactor) delete fMuonScaleFactor;
  if(fChargedHadronScaleFactor) delete fChargedHadronScaleFactor;
}

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
  fCovariance->Calc(fGeometry);
  fCovariance->SetMinHits(fNMinHits);
  // load geometry
  fAcx = fCovariance->AccPnt();

  // import input array
  GetFactory()->EventModel()->Attach(GetString("InputArray", "TrackMerger/tracks"), fInputArray);
  // create output arrays
  GetFactory()->EventModel()->Book(fOutputArray, GetString("OutputArray", "tracks"));
}

//------------------------------------------------------------------------------

void TrackCovariance::Finish()
{
}

//------------------------------------------------------------------------------

void TrackCovariance::Process()
{
  Double_t mass, p, pt, q, ct;
  Double_t dd0, ddz, dphi, dct, dp, dpt, dC;
  //
  // Get cylindrical box for fast track simulation
  //
  Double_t Rin = fGeometry->GetRmin();
  Double_t ZinPos = fGeometry->GetZminPos();
  Double_t ZinNeg = fGeometry->GetZminNeg();

  for(auto &candidate : *fInputArray) //TODO: ensure const-qualification of consumers
  {

    // converting to meters
    auto *particle = static_cast<Candidate *>(candidate.GetCandidates()->At(0));

    // converting to meters
    const auto &candidatePosition = particle->Position * 1e-03;
    const auto &candidateMomentum = particle->Momentum;

    double candidatePositionVect[3], candidateMomentumVect[3];
    candidatePosition.Vect().GetCoordinates(candidatePositionVect);
    candidateMomentum.Vect().GetCoordinates(candidateMomentumVect);

    Bool_t inside = TrkUtil::IsInside(candidatePositionVect, Rin, ZinNeg, ZinPos); // Check if in inner box
    Bool_t Accept = kTRUE;
    if(inside)
      Accept = fCovariance->IsAccepted(candidateMomentumVect);
    else
      Accept = fCovariance->IsAccepted(candidatePositionVect, candidateMomentumVect, fGeometry);
    if(!Accept) continue;

    mass = candidateMomentum.M();

    //
    // Standard implementation with grid
    //ObsTrk track(candidatePosition.Vect(), candidateMomentum.Vect(), candidate.Charge, fCovariance, fGeometry);
    //
    // Try Kalman without any grid
    ObsTrk track(candidatePositionVect, candidateMomentumVect, candidate.Charge, mass, fGeometry);

    // apply rescaling factors to resolution
    if(TMath::Abs(candidate.PID) == 11)
    {
      track.SetScale(fElectronScaleFactor->Eval(candidateMomentum.Pt(), candidateMomentum.Eta(), candidateMomentum.Phi(), candidateMomentum.E(), &candidate));
    }
    else if(TMath::Abs(candidate.PID) == 13)
    {
      track.SetScale(fMuonScaleFactor->Eval(candidateMomentum.Pt(), candidateMomentum.Eta(), candidateMomentum.Phi(), candidateMomentum.E(), &candidate));
    }
    else
    {
      track.SetScale(fChargedHadronScaleFactor->Eval(candidateMomentum.Pt(), candidateMomentum.Eta(), candidateMomentum.Phi(), candidateMomentum.E(), &candidate));
    }

    auto *new_candidate = static_cast<Candidate *>(candidate.Clone());

    const auto track_obs_p = track.GetObsP();
    new_candidate->Momentum.SetCoordinates(track_obs_p.X(), track_obs_p.Y(), track_obs_p.Z(), mass);

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
    new_candidate->P = std::sqrt(track.GetObsP().Mag2());
    new_candidate->PT = pt;
    new_candidate->Charge = q;

    dd0 = TMath::Sqrt(track.GetCov()(0, 0)) * 1e03;
    ddz = TMath::Sqrt(track.GetCov()(3, 3)) * 1e03;
    dphi = TMath::Sqrt(track.GetCov()(1, 1));
    dct = TMath::Sqrt(track.GetCov()(4, 4));
    dpt = 2 * TMath::Sqrt(track.GetCov()(2, 2)) * pt * pt / (0.2998 * fBz);
    dp = TMath::Sqrt((1. + ct * ct) * dpt * dpt + 4 * pt * pt * ct * ct * dct * dct / (1. + ct * ct) / (1. + ct * ct));
    dC = TMath::Sqrt(track.GetCov()(2, 2)) * 1e-03;

    new_candidate->ErrorD0 = dd0;
    new_candidate->ErrorDZ = ddz;
    new_candidate->ErrorP = dp;
    new_candidate->ErrorC = dC;
    new_candidate->ErrorCtgTheta = dct;
    new_candidate->ErrorPhi = dphi;
    new_candidate->ErrorPT = dpt;
    //new_candidate->TrackResolution = dpt / pt;
    new_candidate->TrackResolution = dp / p;

    new_candidate->AddCandidate(&candidate); // mother particle

    fOutputArray->emplace_back(*new_candidate);
  }
}

//------------------------------------------------------------------------------
