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

/** \class ClusterCounting
 *
 *  Counts ionisation clusters of energy loss in drift chambers
 *
 *  \authors F. Bedeschi - INFN Pisa
*            P. Demin - UCLouvain, Louvain-la-Neuve
 *           M. Selvaggi - CERN
 *
 *
 */

#include "TrackCovariance/TrkUtil.h"
#include "classes/DelphesClasses.h"
#include "modules/ClusterCounting.h"

#include "TMath.h"
#include "TVectorD.h"

#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

ClusterCounting::ClusterCounting() :
  fTrackUtil(0)
{
  fTrackUtil = new TrkUtil();
}

//------------------------------------------------------------------------------

ClusterCounting::~ClusterCounting()
{
  if(fTrackUtil) delete fTrackUtil;
}

//------------------------------------------------------------------------------

void ClusterCounting::Init()
{

  // geometric acceptance
  fRmin = GetDouble("Rmin", 0.);
  fRmax = GetDouble("Rmax", 0.);
  fZmin = GetDouble("Zmin", 0.);
  fZmax = GetDouble("Zmax", 0.);

  // magnetic field
  fBz = GetDouble("Bz", 0.);

  // gas mix option: 0
  // 0:  Helium 90 - Isobutane 10
  // 1:  Helium 100
  // 2:  Argon 50 - Ethane 50
  // 3:  Argon 100
  fGasOption = GetInt("GasOption", 0);

  // initialize drift chamber geometry and gas mix
  fTrackUtil->SetBfield(fBz);
  fTrackUtil->SetDchBoundaries(fRmin, fRmax, fZmin, fZmax);
  fTrackUtil->SetGasMix(fGasOption);

  // import input array(s)
  ImportArray(GetString("InputArray", "TrackMerger/tracks"), fInputArray);
  // create output arrays
  ExportArray(fOutputArray, GetString("OutputArray", "tracks"));
}

//------------------------------------------------------------------------------

void ClusterCounting::Finish()
{
}

//------------------------------------------------------------------------------

void ClusterCounting::Process()
{
  Double_t mass, trackLength, Ncl;

  fOutputArray->clear();
  for(const auto &candidate : *fInputArray)
  {
    // converting to meters
    auto *particle = static_cast<Candidate *>(candidate.GetCandidates().at(0));

    // converting to meters
    const auto &candidatePosition = particle->Position * 1e-03;
    const auto &candidateMomentum = particle->Momentum;

    double candidatePositionVect[3], candidateMomentumVect[3];
    candidatePosition.Vect().GetCoordinates(candidatePositionVect);
    candidateMomentum.Vect().GetCoordinates(candidateMomentumVect);

    TVectorD Par = TrkUtil::XPtoPar(TVector3(candidatePositionVect), TVector3(candidateMomentumVect), candidate.Charge, fBz);
    mass = candidateMomentum.M();

    trackLength = fTrackUtil->TrkLen(Par);

    auto new_candidate = candidate;

    Ncl = 0.;
    if(fTrackUtil->IonClusters(Ncl, mass, Par))
    {
      new_candidate.Nclusters = Ncl;
      new_candidate.dNdx = (trackLength > 0.) ? Ncl / trackLength : -1;
    }

    new_candidate.AddCandidate(const_cast<Candidate *>(&candidate)); // preserve parentage

    fOutputArray->emplace_back(new_candidate);
  }
}

//------------------------------------------------------------------------------
