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

/** \class JetPileUpSubtractor
 *
 *  Subtract pile-up contribution from jets using the fastjet area method
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "modules/JetPileUpSubtractor.h"

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

void JetPileUpSubtractor::Init()
{
  fJetPTMin = GetDouble("JetPTMin", 20.0);

  // import input arrays
  GetFactory()->EventModel()->Attach(GetString("JetInputArray", "FastJetFinder/jets"), fJetInputArray);
  GetFactory()->EventModel()->Attach(GetString("RhoInputArray", "Rho/rho"), fRhoInputArray);

  // create output arrays
  GetFactory()->EventModel()->Book(fOutputArray, GetString("OutputArray", "jets"));
}

//------------------------------------------------------------------------------

void JetPileUpSubtractor::Finish()
{
}

//------------------------------------------------------------------------------

void JetPileUpSubtractor::Process()
{
  Double_t eta = 0.0;
  Double_t rho = 0.0;

  // loop over all input candidates
  for(const auto &candidate : *fJetInputArray)
  {
    auto momentum = candidate.Momentum;
    const auto &area = candidate.Area;
    eta = momentum.Eta();

    // find rho
    rho = 0.0;
    if(fRhoInputArray)
    {
      for(const auto &object : *fRhoInputArray)
      {
        if(eta >= object.Edges[0] && eta < object.Edges[1])
        {
          rho = object.Momentum.Pt();
        }
      }
    }

    // apply pile-up correction
    if(momentum.Pt() <= rho * area.Pt()) continue;

    momentum -= rho * area;

    if(momentum.Pt() <= fJetPTMin) continue;

    auto *new_candidate = static_cast<Candidate *>(candidate.Clone());
    new_candidate->Momentum = momentum;

    fOutputArray->emplace_back(*new_candidate);
  }
}

//------------------------------------------------------------------------------
