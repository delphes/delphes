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

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

class AngularSmearing: public DelphesModule
{
public:
  AngularSmearing() :
    fFormulaEta(std::make_unique<DelphesFormula>()), fFormulaPhi(std::make_unique<DelphesFormula>()) {}

  void Init() override;
  void Process() override;

private:
  const std::unique_ptr<DelphesFormula> fFormulaEta; //!
  const std::unique_ptr<DelphesFormula> fFormulaPhi; //!

  CandidatesCollection fInputArray;
  CandidatesCollection fOutputArray;
};

//------------------------------------------------------------------------------

void AngularSmearing::Init()
{
  // read resolution formula
  fFormulaEta->Compile(GetString("EtaResolutionFormula", "0.0"));
  fFormulaPhi->Compile(GetString("PhiResolutionFormula", "0.0"));

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void AngularSmearing::Process()
{
  fOutputArray->clear();
  Double_t pt, eta, phi, e, m;

  for(const auto &candidate : *fInputArray)
  {
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    eta = candidateMomentum.Eta();
    phi = candidateMomentum.Phi();
    pt = candidateMomentum.Pt();
    e = candidateMomentum.E();
    m = candidateMomentum.M();

    // apply smearing formula for eta,phi
    eta = gRandom->Gaus(eta, fFormulaEta->Eval(pt, eta, phi, e, candidate));
    phi = gRandom->Gaus(phi, fFormulaPhi->Eval(pt, eta, phi, e, candidate));

    if(pt <= 0.0) continue;

    auto *new_candidate = static_cast<Candidate *>(candidate->Clone());
    new_candidate->Momentum.SetPtEtaPhiM(pt, eta, phi, m);
    new_candidate->AddCandidate(candidate);

    fOutputArray->emplace_back(new_candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("AngularSmearing", AngularSmearing);
