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

/** \class IdentificationMap
 *
 *  Converts particles with some PDG code into another particle,
 *  according to parametrized probability.
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "modules/IdentificationMap.h"

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

void IdentificationMap::Init()
{
  TMisIDMap::iterator itEfficiencyMap;
  ExRootConfParam param;
  DelphesFormula *formula;
  Int_t i, size, pdg;

  // read efficiency formulas
  param = GetParam("EfficiencyFormula");
  size = param.GetSize();

  fEfficiencyMap.clear();
  for(i = 0; i < size / 3; ++i)
  {
    formula = new DelphesFormula;
    formula->Compile(param[i * 3 + 2].GetString());
    pdg = param[i * 3].GetInt();
    fEfficiencyMap.insert(make_pair(pdg, make_pair(param[i * 3 + 1].GetInt(), formula)));
  }

  // set default efficiency formula
  itEfficiencyMap = fEfficiencyMap.find(0);
  if(itEfficiencyMap == fEfficiencyMap.end())
  {
    formula = new DelphesFormula;
    formula->Compile("1.0");

    fEfficiencyMap.insert(make_pair(0, make_pair(0, formula)));
  }

  // import input array
  GetFactory()->EventModel()->Attach(GetString("InputArray", "ParticlePropagator/stableParticles"), fInputArray);
  // create output array
  ExportArray(fOutputArray, GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void IdentificationMap::Finish()
{
  for(auto itEfficiencyMap = fEfficiencyMap.begin(); itEfficiencyMap != fEfficiencyMap.end(); ++itEfficiencyMap)
    if(auto *formula = (itEfficiencyMap->second).second; formula) delete formula;
}

//------------------------------------------------------------------------------

void IdentificationMap::Process()
{
  Double_t pt, eta, phi, e;
  TMisIDMap::iterator itEfficiencyMap;
  pair<TMisIDMap::iterator, TMisIDMap::iterator> range;
  DelphesFormula *formula;
  Int_t pdgCodeIn, pdgCodeOut, charge;

  Double_t p, r, total;

  for(const auto &candidate : *fInputArray)
  {
    const auto &candidatePosition = candidate.Position;
    const auto &candidateMomentum = candidate.Momentum;
    eta = candidatePosition.Eta();
    phi = candidatePosition.Phi();
    pt = candidateMomentum.Pt();
    e = candidateMomentum.E();

    pdgCodeIn = candidate.PID;
    charge = candidate.Charge;

    // first check that PID of this particle is specified in the map
    // otherwise, look for PID = 0

    itEfficiencyMap = fEfficiencyMap.find(pdgCodeIn);

    range = fEfficiencyMap.equal_range(pdgCodeIn);
    if(range.first == range.second) range = fEfficiencyMap.equal_range(-pdgCodeIn);
    if(range.first == range.second) range = fEfficiencyMap.equal_range(0);

    r = gRandom->Uniform();
    total = 0.0;

    // loop over sub-map for this PID
    for(TMisIDMap::iterator it = range.first; it != range.second; ++it)
    {
      formula = (it->second).second;
      pdgCodeOut = (it->second).first;

      p = formula->Eval(pt, eta, phi, e);

      if(total <= r && r < total + p)
      {
        // change PID of particle
        auto new_candidate = candidate;
        if(pdgCodeOut != 0) new_candidate.PID = charge * pdgCodeOut;
        fOutputArray->emplace_back(new_candidate);
        break;
      }

      total += p;
    }
  }
}

//------------------------------------------------------------------------------
