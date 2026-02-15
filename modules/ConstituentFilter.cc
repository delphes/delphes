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

/** \class ConstituentFilter
 *
 *  Drops all input objects that are not constituents of any jet.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/ConstituentFilter.h"

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

void ConstituentFilter::Init()
{
  ExRootConfParam param;

  fJetPTMin = GetDouble("JetPTMin", 0.0);

  // import input array(s)

  param = GetParam("JetInputArray");
  for(Long_t i = 0; i < param.GetSize(); ++i)
    ImportArray(param[i].GetString(), fInputList.emplace_back());

  param = GetParam("ConstituentInputArray");
  for(Long_t i = 0; i < param.GetSize() / 2; ++i)
  {
    auto &[input_collection, output_collection] = fInputMap.emplace_back();
    ImportArray(param[i * 2].GetString(), input_collection);
    ExportArray(output_collection, param[i * 2 + 1].GetString());
  }
}

//------------------------------------------------------------------------------

void ConstituentFilter::Finish()
{
}

//------------------------------------------------------------------------------

void ConstituentFilter::Process()
{
  for(const auto &[input_collection, output_collection] : fInputMap)
    output_collection->clear();

  // loop over all jet input arrays
  for(const auto &input_collection : fInputList)
  {
    // loop over all jets
    for(auto &jet : *input_collection) //TODO: ensure cons-qualification
    {
      if(jet.Momentum.Pt() <= fJetPTMin) continue;

      // loop over all constituents
      for(const auto &constituent : jet.GetCandidates())
      {
        // set the IsConstituent flag
        constituent->IsConstituent = 1;
      }
    }
  }

  // loop over all constituent input arrays
  for(const auto &[input_collection, output_collection] : fInputMap)
  {
    // loop over all constituents
    for(const auto &constituent : *input_collection)
      if(constituent.IsConstituent) // check the IsConstituent flag
        output_collection->emplace_back(constituent);
  }
}

//------------------------------------------------------------------------------
