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

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm> 
#include <stdexcept>
#include <iostream>
#include <sstream>

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
  // import arrays with output from other modules

  ExRootConfParam param = GetParam("InputArray");
  Long_t i, size;
  const TObjArray *array;
  TIterator *iterator;

  fInputMap.clear();

  size = param.GetSize();
  for(i = 0; i < size/2; ++i)
  {
    array = ImportArray(param[i*2].GetString());
    iterator = array->MakeIterator();

    fInputMap.push_back(make_pair(iterator, ExportArray(param[i*2 + 1].GetString())));
  }
}

//------------------------------------------------------------------------------

void UniqueObjectFinder::Finish()
{
  vector< pair< TIterator *, TObjArray * > >::iterator itInputMap;
  TIterator *iterator;

  for(itInputMap = fInputMap.begin(); itInputMap != fInputMap.end(); ++itInputMap)
  {
    iterator = itInputMap->first;

    if(iterator) delete iterator;
  }
}

//------------------------------------------------------------------------------

void UniqueObjectFinder::Process()
{
  Candidate *candidate;
  vector< pair< TIterator *, TObjArray * > >::iterator itInputMap;
  TIterator *iterator;
  TObjArray *array;

  // loop over all input arrays
  for(itInputMap = fInputMap.begin(); itInputMap != fInputMap.end(); ++itInputMap)
  {
    iterator = itInputMap->first;
    array = itInputMap->second;

    // loop over all candidates
    iterator->Reset();
    while((candidate = static_cast<Candidate*>(iterator->Next())))
    {
      if(Unique(candidate, itInputMap))
      {
        array->Add(candidate);
      }
    }
  }
}

//------------------------------------------------------------------------------

Bool_t UniqueObjectFinder::Unique(Candidate *candidate, vector< pair< TIterator *, TObjArray * > >::iterator itInputMap)
{
  Candidate *previousCandidate;
  vector< pair< TIterator *, TObjArray * > >::iterator previousItInputMap;
  TObjArray *array;

  // loop over previous arrays
  for(previousItInputMap = fInputMap.begin(); previousItInputMap != itInputMap; ++previousItInputMap)
  {
    array = previousItInputMap->second;
    TIter iterator(array);

    // loop over all candidates
    iterator.Reset();
    while((previousCandidate = static_cast<Candidate*>(iterator.Next())))
    {
      if(candidate->Overlaps(previousCandidate))
      {
        return kFALSE;
      }
    }
  }

  return kTRUE;
}

//------------------------------------------------------------------------------
