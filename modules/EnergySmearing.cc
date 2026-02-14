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

/** \class EnergySmearing
 *
 *  Performs energy resolution smearing.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/EnergySmearing.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

EnergySmearing::EnergySmearing() :
  fFormula(0)
{
  fFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

EnergySmearing::~EnergySmearing()
{
  if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void EnergySmearing::Init()
{
  // read resolution formula

  fFormula->Compile(GetString("ResolutionFormula", "0.0"));

  // import input array(s)
  ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"), fInputArray);
  // create output arrays
  ExportArray(fOutputArray, GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void EnergySmearing::Finish()
{
}

//------------------------------------------------------------------------------

void EnergySmearing::Process()
{
  Double_t pt, energy, eta, phi, m;

  for(const auto &candidate : *fInputArray)
  {
    const auto &candidatePosition = candidate.Position;
    const auto &candidateMomentum = candidate.Momentum;

    pt = candidatePosition.Pt();
    eta = candidatePosition.Eta();
    phi = candidatePosition.Phi();
    energy = candidateMomentum.E();
    m = candidateMomentum.M();

    // apply smearing formula
    energy = gRandom->Gaus(energy, fFormula->Eval(pt, eta, phi, energy));

    if(energy <= 0.0) continue;

    auto new_candidate = candidate;
    eta = candidateMomentum.Eta();
    phi = candidateMomentum.Phi();
    pt = (energy > m) ? TMath::Sqrt(energy * energy - m * m) / TMath::CosH(eta) : 0;
    new_candidate.Momentum = ROOT::Math::PtEtaPhiEVector(pt, eta, phi, energy);
    new_candidate.TrackResolution = fFormula->Eval(pt, eta, phi, energy) / candidateMomentum.E();
    new_candidate.AddCandidate(const_cast<Candidate *>(&candidate)); //TODO: ensure const-qualification
    fOutputArray->emplace_back(new_candidate);
  }
}

//------------------------------------------------------------------------------
