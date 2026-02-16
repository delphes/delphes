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

/** \class AngularSmearing
 *
 *  Performs transverse angular resolution smearing.
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "modules/AngularSmearing.h"

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

AngularSmearing::AngularSmearing() :
  fFormulaEta(0), fFormulaPhi(0)
{
  fFormulaEta = new DelphesFormula;
  fFormulaPhi = new DelphesFormula;
}

//------------------------------------------------------------------------------

AngularSmearing::~AngularSmearing()
{
  if(fFormulaEta) delete fFormulaEta;
  if(fFormulaPhi) delete fFormulaPhi;
}

//------------------------------------------------------------------------------

void AngularSmearing::Init()
{
  // read resolution formula

  fFormulaEta->Compile(GetString("EtaResolutionFormula", "0.0"));
  fFormulaPhi->Compile(GetString("PhiResolutionFormula", "0.0"));

  // import input array
  fInputArray = ImportArray<std::vector<Candidate> >(GetString("InputArray", "ParticlePropagator/stableParticles"));

  // create output array
  fOutputArray = ExportArray<std::vector<Candidate> >(GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void AngularSmearing::Finish()
{
}

//------------------------------------------------------------------------------

void AngularSmearing::Process()
{
  Double_t pt, eta, phi, e, m;
  fOutputArray->clear();
  for(const auto &candidate : *fInputArray)
  {
    const auto &candidateMomentum = candidate.Momentum;
    eta = candidateMomentum.Eta();
    phi = candidateMomentum.Phi();
    pt = candidateMomentum.Pt();
    e = candidateMomentum.E();
    m = candidateMomentum.M();

    // apply smearing formula for eta,phi
    eta = gRandom->Gaus(eta, fFormulaEta->Eval(pt, eta, phi, e, const_cast<Candidate *>(&candidate))); //TODO: const-qualified version?
    phi = gRandom->Gaus(phi, fFormulaPhi->Eval(pt, eta, phi, e, const_cast<Candidate *>(&candidate)));

    if(pt <= 0.0) continue;

    auto new_candidate = candidate;
    new_candidate.Momentum = ROOT::Math::PtEtaPhiMVector(pt, eta, phi, m);
    new_candidate.AddCandidate(&candidate); // preserve parentage

    fOutputArray->emplace_back(new_candidate);
  }
}

//------------------------------------------------------------------------------
