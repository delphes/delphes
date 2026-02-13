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

//------------------------------------------------------------------------------

/** \class TaggingParticlesSkimmer
 *
 *  Filters particle collection by only keeping gen particles useful for b and tau tagging.
 *  tau particles are replaced by their "visible decay".
 *
 *  \author M. Selvaggi
 *
 */

#include "modules/TaggingParticlesSkimmer.h"
#include "modules/TauTagging.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootSTLVectorFilter.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TLorentzVector.h"
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

TaggingParticlesSkimmer::TaggingParticlesSkimmer() :
  fClassifier(0), fFilter(0)
{
}

//------------------------------------------------------------------------------

TaggingParticlesSkimmer::~TaggingParticlesSkimmer()
{
}

//------------------------------------------------------------------------------

void TaggingParticlesSkimmer::Init()
{

  fPTMin = GetDouble("PTMin", 15.0);
  fEtaMax = GetDouble("EtaMax", 2.5);

  // import input arrays
  GetFactory()->EventModel()->Attach(GetString("InputArray", "Delphes/partons"), fPartonInputArray);
  GetFactory()->EventModel()->Attach(GetString("InputArray", "Delphes/allParticles"), fParticleInputArray);

  // create output array
  GetFactory()->EventModel()->Book(fOutputArray, GetString("OutputArray", "taggingParticles"));

  fClassifier = new TauTaggingPartonClassifier(*fParticleInputArray);
  fClassifier->fPTMin = GetDouble("PTMin", 15.0);
  fClassifier->fEtaMax = GetDouble("EtaMax", 2.5);

  fFilter = new ExRootSTLVectorFilter(*fPartonInputArray);
}

//------------------------------------------------------------------------------

void TaggingParticlesSkimmer::Finish()
{
  if(fFilter) delete fFilter;
  if(fClassifier) delete fClassifier;
}

//------------------------------------------------------------------------------

void TaggingParticlesSkimmer::Process()
{
  TLorentzVector tauMomentum;
  Double_t pt, eta;
  Int_t pdgCode;

  // first select hadronic taus and replace them by visible part
  fFilter->Reset();
  const auto tauArray = fFilter->GetSubArray(fClassifier, 0);

  if(tauArray.empty()) return;

  // loop over all input taus
  for(const auto &tau : tauArray)
  {
    if(tau.D1 < 0) continue;

    if(tau.D1 >= static_cast<int>(fParticleInputArray->size()) || tau.D2 >= static_cast<int>(fParticleInputArray->size()))
    {
      throw runtime_error("tau's daughter index is greater than the ParticleInputArray size");
    }

    tauMomentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);

    for(int i = tau.D1; i <= tau.D2; ++i)
    {
      const auto &daughter = fParticleInputArray->at(i);
      if(TMath::Abs(daughter.PID) == 16) continue;
      tauMomentum += daughter.Momentum;
    }

    auto *new_candidate = static_cast<Candidate *>(tau.Clone());
    new_candidate->Momentum = tauMomentum;

    fOutputArray->emplace_back(*new_candidate);
  }

  // then add all other partons (except tau's to avoid double counting)
  for(const auto &candidate : *fPartonInputArray)
  {
    pdgCode = TMath::Abs(candidate.PID);
    if(pdgCode == 15) continue;

    pt = candidate.Momentum.Pt();
    if(pt < fPTMin) continue;

    eta = TMath::Abs(candidate.Momentum.Eta());
    if(eta > fEtaMax) continue;

    fOutputArray->emplace_back(candidate);
  }
}
