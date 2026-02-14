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

/** \class Merger
 *
 *  Merges multiple input arrays into one output array
 *  and sums transverse momenta of all input objects.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/Merger.h"

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

void Merger::Init()
{
  // import arrays with output from other modules

  ExRootConfParam param = GetParam("InputArray");

  for(Long_t i = 0; i < param.GetSize(); ++i)
    GetFactory()->EventModel()->Attach(param[i].GetString(), fInputList.emplace_back());

  // create output arrays
  GetFactory()->EventModel()->Book(fOutputArray, GetString("OutputArray", "candidates"));
  GetFactory()->EventModel()->Book(fMomentumOutputArray, GetString("MomentumOutputArray", "momentum"));
  GetFactory()->EventModel()->Book(fEnergyOutputArray, GetString("EnergyOutputArray", "energy"));
}

//------------------------------------------------------------------------------

void Merger::Finish()
{
}

//------------------------------------------------------------------------------

void Merger::Process()
{
  Double_t sumPT, sumE;

  DelphesFactory *factory = GetFactory();

  ROOT::Math::PxPyPzEVector momentum;
  sumPT = 0;
  sumE = 0;

  // loop over all input arrays
  for(const auto &input_collection : fInputList)
  {
    // loop over all candidates
    for(const auto &candidate : *input_collection)
    {
      const auto &candidateMomentum = candidate.Momentum;

      momentum += candidateMomentum;
      sumPT += candidateMomentum.Pt();
      sumE += candidateMomentum.E();

      fOutputArray->emplace_back(candidate);
    }
  }

  auto *momentum_candidate = factory->NewCandidate();
  momentum_candidate->Position.SetXYZT(0.0, 0.0, 0.0, 0.0);
  momentum_candidate->Momentum = momentum;
  fMomentumOutputArray->emplace_back(*momentum_candidate);

  auto *energy_candidate = factory->NewCandidate();
  energy_candidate->Position.SetXYZT(0.0, 0.0, 0.0, 0.0);
  energy_candidate->Momentum = ROOT::Math::PtEtaPhiEVector(sumPT, 0.0, 0.0, sumE);
  fEnergyOutputArray->emplace_back(*energy_candidate);
}

//------------------------------------------------------------------------------
