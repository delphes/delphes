/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
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

/** \class BeamSpotSmearing
 *
 *  Applies per-event Gaussian smearing to the beam spot and the bunch time spread.
 *
 *  \author Michele Selvaggi - CERN
 *
 */

#include "modules/BeamSpotSmearing.h"

#include "classes/DelphesClasses.h"

#include "TIterator.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TRandom3.h"

using namespace std;

//------------------------------------------------------------------------------

BeamSpotSmearing::BeamSpotSmearing() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

BeamSpotSmearing::~BeamSpotSmearing()
{
}

//------------------------------------------------------------------------------

void BeamSpotSmearing::Init()
{
  const Double_t c_light = 2.99792458E8;

  // beam spot size in x, y, z given in [m], bunch time spread in [s];
  fSigmaX = GetDouble("SigmaX", 0.0) * 1.0E3;
  fSigmaY = GetDouble("SigmaY", 0.0) * 1.0E3;
  fSigmaZ = GetDouble("SigmaZ", 0.0) * 1.0E3;
  fSigmaT = GetDouble("SigmaT", 0.0) * c_light * 1.0E3;

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));
  fItInputArray = fInputArray->MakeIterator();
}

//------------------------------------------------------------------------------

void BeamSpotSmearing::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void BeamSpotSmearing::Process()
{
  Candidate *candidate;
  Double_t dx, dy, dz, dt;

  // do nothing unless at least one spread is configured
  if(fSigmaX <= 0.0 && fSigmaY <= 0.0 && fSigmaZ <= 0.0 && fSigmaT <= 0.0) return;

  // one rigid shift per event; draws are skipped for vanishing spreads so
  // that they consume no random numbers
  dx = fSigmaX > 0.0 ? gRandom->Gaus(0.0, fSigmaX) : 0.0;
  dy = fSigmaY > 0.0 ? gRandom->Gaus(0.0, fSigmaY) : 0.0;
  dz = fSigmaZ > 0.0 ? gRandom->Gaus(0.0, fSigmaZ) : 0.0;
  dt = fSigmaT > 0.0 ? gRandom->Gaus(0.0, fSigmaT) : 0.0;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    TLorentzVector &position = candidate->Position;
    position.SetXYZT(position.X() + dx, position.Y() + dy, position.Z() + dz, position.T() + dt);

    TLorentzVector &decayPosition = candidate->DecayPosition;
    decayPosition.SetXYZT(decayPosition.X() + dx, decayPosition.Y() + dy, decayPosition.Z() + dz, decayPosition.T() + dt);
  }
}

//------------------------------------------------------------------------------
