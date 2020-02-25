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
 *  \author P. Demin - UCLouvain, Louvain-la-Neuve
 *
 */

//FIXME add reference to Bedeschi-code
//FIXME add smearing of impact parameter as well 


#include "modules/TrackCovariance.h"

#include "classes/DelphesClasses.h"

#include "TrackCovariance/SolGeom.h"
#include "TrackCovariance/SolGridCov.h"
#include "TrackCovariance/ObsTrk.h"

#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"

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
  Double_t energy;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;

    energy = candidateMomentum.E();

    ObsTrk track(candidatePosition.Vect(), candidateMomentum.Vect(), candidate->Charge, fBz, fCovariance);

    mother = candidate;
    candidate = static_cast<Candidate *>(candidate->Clone());
    candidate->Momentum.SetVect(track.GetObsP());
    candidate->Momentum.SetE(energy);

    candidate->AddCandidate(mother);

    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
