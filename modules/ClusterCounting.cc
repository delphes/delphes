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

#include "modules/ClusterCounting.h"
#include "classes/DelphesClasses.h"
#include "TrackCovariance/TrkUtil.h"

#include "TLorentzVector.h"
#include "TVectorD.h"
#include "TMath.h"
#include "TObjArray.h"

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
  fRmin  = GetDouble("Rmin", 0.);
  fRmax  = GetDouble("Rmax", 0.);
  fZmin  = GetDouble("Zmin", 0.);
  fZmax  = GetDouble("Zmax", 0.);

  // magnetic field
  fBz    = GetDouble("Bz", 0.);

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

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "TrackMerger/tracks"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
}

//------------------------------------------------------------------------------

void ClusterCounting::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void ClusterCounting::Process()
{
  Candidate *candidate, *mother, *particle;
  Double_t mass, trackLength, Ncl;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {

    // converting to meters
    particle = static_cast<Candidate *>(candidate->GetCandidates()->At(0));

    // converting to meters
    const TLorentzVector &candidatePosition = particle->Position*1e-03;
    const TLorentzVector &candidateMomentum = particle->Momentum;

		TVectorD Par = TrkUtil::XPtoPar(candidatePosition.Vect(), candidateMomentum.Vect(), candidate->Charge, fBz);
    mass = candidateMomentum.M();

    trackLength = fTrackUtil->TrkLen(Par);

    mother    = candidate;
    candidate = static_cast<Candidate*>(candidate->Clone());

    Ncl = -999;


    // computation of Nclusters is not supported for electrons
    /*
    if (TMath::Abs(particle->PID) == 11)
    {
      candidate->Nclusters = Ncl;
      candidate->dNdx = -999;
    }
    */
    if (fTrackUtil->IonClusters(Ncl, mass, Par))
    {
      candidate->Nclusters = Ncl;
      candidate->dNdx = (trackLength > 0.) ? Ncl/trackLength : -1;
    }

    candidate->AddCandidate(mother);

    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
