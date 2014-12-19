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

ConstituentFilter::ConstituentFilter()
{
}

//------------------------------------------------------------------------------

ConstituentFilter::~ConstituentFilter()
{
}

//------------------------------------------------------------------------------

void ConstituentFilter::Init()
{
  ExRootConfParam param;
  Long_t i, size;
  const TObjArray *array;
  TIterator *iterator;

  fJetPTMin = GetDouble("JetPTMin", 0.0);

  // import input array(s)

  param = GetParam("JetInputArray");
  size = param.GetSize();
  for(i = 0; i < size; ++i)
  {
    array = ImportArray(param[i].GetString());
    iterator = array->MakeIterator();

    fInputList.push_back(iterator);
  }

  param = GetParam("ConstituentInputArray");
  size = param.GetSize();
  for(i = 0; i < size/2; ++i)
  {
    array = ImportArray(param[i*2].GetString());
    iterator = array->MakeIterator();

    fInputMap[iterator] = ExportArray(param[i*2 + 1].GetString());
  }
}

//------------------------------------------------------------------------------

void ConstituentFilter::Finish()
{
  map< TIterator *, TObjArray * >::iterator itInputMap;
  vector< TIterator * >::iterator itInputList;
  TIterator *iterator;

  for(itInputList = fInputList.begin(); itInputList != fInputList.end(); ++itInputList)
  {
    iterator = *itInputList;
    if(iterator) delete iterator;
  }

  for(itInputMap = fInputMap.begin(); itInputMap != fInputMap.end(); ++itInputMap)
  {
    iterator = itInputMap->first;
    if(iterator) delete iterator;
  }
}

//------------------------------------------------------------------------------

void ConstituentFilter::Process()
{
  Candidate *jet, *constituent;
  map< TIterator *, TObjArray * >::iterator itInputMap;
  vector< TIterator * >::iterator itInputList;
  TIterator *iterator;
  TObjArray *array;

  // loop over all jet input arrays
  for(itInputList = fInputList.begin(); itInputList != fInputList.end(); ++itInputList)
  {
    iterator = *itInputList;

    // loop over all jets
    iterator->Reset();
    while((jet = static_cast<Candidate*>(iterator->Next())))
    {
      TIter itConstituents(jet->GetCandidates());

      if(jet->Momentum.Pt() <= fJetPTMin) continue;

      // loop over all constituents
      while((constituent = static_cast<Candidate*>(itConstituents.Next())))
      {
        // set the IsConstituent flag
        constituent->IsConstituent = 1;
      }
    }
  }

  // loop over all constituent input arrays
  for(itInputMap = fInputMap.begin(); itInputMap != fInputMap.end(); ++itInputMap)
  {
    iterator = itInputMap->first;
    array = itInputMap->second;

    // loop over all constituents
    iterator->Reset();
    while((constituent = static_cast<Candidate*>(iterator->Next())))
    {
      // check the IsConstituent flag
      if(constituent->IsConstituent)
      {
        array->Add(constituent);
      }
    }
  }
}

//------------------------------------------------------------------------------
