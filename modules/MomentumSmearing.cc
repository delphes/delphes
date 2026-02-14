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

/** \class MomentumSmearing
 *
 *  Performs transverse momentum resolution smearing.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/MomentumSmearing.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

MomentumSmearing::MomentumSmearing() :
  fFormula(0)
{
  fFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

MomentumSmearing::~MomentumSmearing()
{
  if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void MomentumSmearing::Init()
{
  // read resolution formula

  fFormula->Compile(GetString("ResolutionFormula", "0.0"));

  // import input array
  GetFactory()->EventModel()->Attach(GetString("InputArray", "ParticlePropagator/stableParticles"), fInputArray);

  // switch to compute momentum smearing based on momentum vector eta, phi
  fUseMomentumVector = GetBool("UseMomentumVector", false);

  // create output array
  GetFactory()->EventModel()->Book(fOutputArray, GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void MomentumSmearing::Finish()
{
}

//------------------------------------------------------------------------------

void MomentumSmearing::Process()
{
  Double_t pt, eta, phi, e, m, res;

  for(const auto &candidate : *fInputArray)
  {
    const auto &candidatePosition = candidate.Position;
    const auto &candidateMomentum = candidate.Momentum;
    eta = candidatePosition.Eta();
    phi = candidatePosition.Phi();

    if(fUseMomentumVector)
    {
      eta = candidateMomentum.Eta();
      phi = candidateMomentum.Phi();
    }

    pt = candidateMomentum.Pt();
    e = candidateMomentum.E();
    m = candidateMomentum.M();
    res = fFormula->Eval(pt, eta, phi, e, const_cast<Candidate *>(&candidate)); //TODO: check if we can use const qualifier

    // apply smearing formula
    //pt = gRandom->Gaus(pt, fFormula->Eval(pt, eta, phi, e) * pt);

    res = (res > 1.0) ? 1.0 : res;

    pt = LogNormal(pt, res * pt);

    //if(pt <= 0.0) continue;

    auto *new_candidate = static_cast<Candidate *>(candidate.Clone());
    eta = candidateMomentum.Eta();
    phi = candidateMomentum.Phi();
    new_candidate->Momentum = ROOT::Math::PtEtaPhiMVector(pt, eta, phi, m);
    //new_candidate->TrackResolution = fFormula->Eval(pt, eta, phi, e);
    new_candidate->TrackResolution = res;
    new_candidate->AddCandidate(const_cast<Candidate *>(&candidate)); // ensure parentage

    fOutputArray->emplace_back(*new_candidate);
  }
}
//----------------------------------------------------------------

Double_t MomentumSmearing::LogNormal(Double_t mean, Double_t sigma)
{
  Double_t a, b;

  if(mean > 0.0)
  {
    b = TMath::Sqrt(TMath::Log((1.0 + (sigma * sigma) / (mean * mean))));
    a = TMath::Log(mean) - 0.5 * b * b;

    return TMath::Exp(a + b * gRandom->Gaus(0.0, 1.0));
  }
  else
  {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
