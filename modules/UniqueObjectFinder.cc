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

/** \class UniqueObjectFinder
 *
 *  Finds uniquely identified photons, electrons and jets.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/UniqueObjectFinder.h"

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

//TODO: large refactoring of this code with change of memory management. A thorough validation is required!

//------------------------------------------------------------------------------

void UniqueObjectFinder::Init()
{
  // use GetUniqueID algorithm to find unique objects (faster than the default Overlaps method)
  fUseUniqueID = GetBool("UseUniqueID", false);

  // import arrays with output from other modules

  ExRootConfParam param = GetParam("InputArray");
  Long_t i, size;

  fInputMap.clear();

  size = param.GetSize();
  for(i = 0; i < size / 2; ++i)
  {
    auto &[input_collection, output_collection] = fInputMap.emplace_back();
    ImportArray(param[i * 2].GetString(), input_collection);
    ExportArray(output_collection, param[i * 2 + 1].GetString());
  }
}

//------------------------------------------------------------------------------

void UniqueObjectFinder::Finish()
{
}

//------------------------------------------------------------------------------

void UniqueObjectFinder::Process()
{
  // loop over all input arrays
  for(const auto &[input_collection, output_collection] : fInputMap)
  {
    output_collection->clear();
    // loop over all candidates
    for(const auto &candidate : *input_collection)
    {
      if(Unique(&candidate, fInputMap.begin()))
      {
        output_collection->emplace_back(candidate);
      }
    }
  }
}

//------------------------------------------------------------------------------

Bool_t UniqueObjectFinder::Unique(const Candidate *candidate, InputMap::iterator itInputMap)
{
  // loop over previous arrays
  for(auto &[input_collection, output_collection] : fInputMap)
  {
    // loop over all candidates
    for(const auto &other_candidate : *output_collection)
    {
      if(fUseUniqueID)
      {
        if(candidate->GetUniqueID() == other_candidate.GetUniqueID())
        {
          return kFALSE;
        }
      }
      else
      {
        if(candidate->Overlaps(&other_candidate))
        {
          return kFALSE;
        }
      }
    }
  }

  return kTRUE;
}

//------------------------------------------------------------------------------
