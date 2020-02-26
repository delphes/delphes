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

#include "modules/TrackCovariance.h"

#include "classes/DelphesClasses.h"

#include "TrackCovariance/SolGeom.h"
#include "TrackCovariance/SolGridCov.h"
#include "TrackCovariance/ObsTrk.h"

#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"

#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

TrackCovariance::TrackCovariance() :
  fGeometry(0), fCovariance(0), fItInputArray(0)
{
  fGeometry = new SolGeom();
  fCovariance = new SolGridCov();
}

//------------------------------------------------------------------------------

TrackCovariance::~TrackCovariance()
{
  if(fGeometry) delete fGeometry;
  if(fCovariance) delete fCovariance;
}

//------------------------------------------------------------------------------

void TrackCovariance::Init()
{
  fBz = GetDouble("Bz", 0.0);
  fGeometry->Read(GetString("DetectorGeometry", ""));

  fCovariance->Calc(fGeometry);

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "TrackMerger/tracks"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
}

//------------------------------------------------------------------------------

void TrackCovariance::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void TrackCovariance::Process()
{
  Candidate *candidate, *mother;
  Double_t mass, p, pt, q, ct;
  Double_t dd0, ddz, dphi, dct, dp, dpt;
  

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    const TLorentzVector &candidatePosition = candidate->InitialPosition;
    const TLorentzVector &candidateMomentum = candidate->Momentum;

    mass = candidateMomentum.M();

    ObsTrk track(candidatePosition.Vect(), candidateMomentum.Vect(), candidate->Charge, fBz, fCovariance);

    mother = candidate;
    candidate = static_cast<Candidate *>(candidate->Clone());

    candidate->Momentum.SetVectM(track.GetObsP(), mass);
    candidate->InitialPosition.SetXYZT(track.GetObsX().X(),track.GetObsX().Y(),track.GetObsX().Z(),candidatePosition.T());

    pt = candidate->Momentum.Pt();
    p  = candidate->Momentum.P();
    q  = track.GetObsQ();
    ct = track.GetObsPar()[4];

    candidate->D0 = track.GetObsPar()[0];
    candidate->DZ = track.GetObsPar()[3];
    candidate->P  = track.GetObsP().Mag();
    candidate->CtgTheta = track.GetObsPar()[4];
    candidate->Phi = track.GetObsPar()[1];

    candidate->PT = pt;
    candidate->Charge = q;

    dd0       = TMath::Sqrt(track.GetCov()(0, 0)); 
    ddz       = TMath::Sqrt(track.GetCov()(3, 3)); 
    dphi      = TMath::Sqrt(track.GetCov()(1, 1)); 
    dct       = TMath::Sqrt(track.GetCov()(4, 4)); 
    dpt       = 2 * TMath::Sqrt( track.GetCov()(2, 2))*pt*pt / (0.2998*fBz);
    dp        = TMath::Sqrt((1.+ct*ct)*dpt*dpt + 4*pt*pt*ct*ct*dct*dct/(1.+ct*ct)/(1.+ct*ct));

    candidate->ErrorD0 = dd0;
    candidate->ErrorDZ = ddz;
    candidate->ErrorP = dp;
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
