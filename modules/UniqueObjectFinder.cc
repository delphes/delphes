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
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

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
#include <utility>

using namespace std;

//------------------------------------------------------------------------------

UniqueObjectFinder::UniqueObjectFinder()
{
}

//------------------------------------------------------------------------------

UniqueObjectFinder::~UniqueObjectFinder()
{
}

//------------------------------------------------------------------------------

void UniqueObjectFinder::Init()
{
  // use GetUniqueID algorithm to find unique objects (faster than the default Overlaps method)
  fUseUniqueID = GetBool("UseUniqueID", false);

  // import arrays with output from other modules

  ExRootConfParam param = GetParam("InputArray");
  Long_t i, size;
  const TObjArray *array;
  TEntryStruct entry;

  size = param.GetSize();

  fInputList.clear();
  for(i = 0; i < size / 2; ++i)
  {
    array = ImportArray(param[i * 2].GetString());
    entry.iterator.reset(array->MakeIterator());
    entry.array = ExportArray(param[i * 2 + 1].GetString());
    fInputList.push_back(move(entry));
  }
}

//------------------------------------------------------------------------------

void UniqueObjectFinder::Finish()
{
}

//------------------------------------------------------------------------------

void UniqueObjectFinder::Process()
{
  Candidate *candidate;
  vector<TEntryStruct>::iterator itInputList;
  TIterator *iterator;
  TObjArray *array;

  // loop over all input arrays
  for(itInputList = fInputList.begin(); itInputList != fInputList.end(); ++itInputList)
  {
    iterator = itInputList->iterator.get();
    array = itInputList->array;

    // loop over all candidates
    iterator->Reset();
    while((candidate = static_cast<Candidate *>(iterator->Next())))
    {
      if(Unique(candidate, itInputList))
      {
        array->Add(candidate);
      }
    }
  }
}

//------------------------------------------------------------------------------

Bool_t UniqueObjectFinder::Unique(Candidate *candidate, vector<TEntryStruct>::iterator itInputList)
{
  Candidate *previousCandidate;
  vector<TEntryStruct>::iterator previousItInputList;
  TObjArray *array;

  // loop over previous arrays
  for(previousItInputList = fInputList.begin(); previousItInputList != itInputList; ++previousItInputList)
  {
    array = previousItInputList->array;
    TIter iterator(array);

    // loop over all candidates
    iterator.Reset();
    while((previousCandidate = static_cast<Candidate *>(iterator.Next())))
    {
      if(fUseUniqueID)
      {
        if(candidate->GetUniqueID() == previousCandidate->GetUniqueID())
        {
          return kFALSE;
        }
      }
      else
      {
        if(candidate->Overlaps(previousCandidate))
        {
          return kFALSE;
        }
      }
    }
  }

  return kTRUE;
}

//------------------------------------------------------------------------------
