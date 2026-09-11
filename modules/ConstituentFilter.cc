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
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

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
  TObjArray *array;
  TEntryStruct entry;

  fJetPTMin = GetDouble("JetPTMin", 0.0);

  // import input array(s)

  param = GetParam("JetInputArray");
  size = param.GetSize();

  fJetList.clear();
  for(i = 0; i < size; ++i)
  {
    array = ImportArray(param[i].GetString());
    entry.iterator.reset(array->MakeIterator());
    entry.array = array;
    fJetList.push_back(move(entry));
  }

  param = GetParam("ConstituentInputArray");
  size = param.GetSize();

  fConstituentList.clear();
  for(i = 0; i < size / 2; ++i)
  {
    array = ImportArray(param[i * 2].GetString());
    entry.iterator.reset(array->MakeIterator());
    entry.array = ExportArray(param[i * 2 + 1].GetString());
    fConstituentList.push_back(move(entry));
  }
}

//------------------------------------------------------------------------------

void ConstituentFilter::Finish()
{
}

//------------------------------------------------------------------------------

void ConstituentFilter::Process()
{
  Candidate *jet, *constituent;
  vector<TEntryStruct>::iterator itList;
  TIterator *iterator;
  TObjArray *array;

  // loop over all jet input arrays
  for(itList = fJetList.begin(); itList != fJetList.end(); ++itList)
  {
    iterator = itList->iterator.get();

    // loop over all jets
    iterator->Reset();
    while((jet = static_cast<Candidate *>(iterator->Next())))
    {
      TIter itConstituents(jet->GetCandidates());

      if(jet->Momentum.Pt() <= fJetPTMin) continue;

      // loop over all constituents
      while((constituent = static_cast<Candidate *>(itConstituents.Next())))
      {
        // set the IsConstituent flag
        constituent->IsConstituent = 1;
      }
    }
  }

  // loop over all constituent input arrays
  for(itList = fConstituentList.begin(); itList != fConstituentList.end(); ++itList)
  {
    iterator = itList->iterator.get();
    array = itList->array;

    // loop over all constituents
    iterator->Reset();
    while((constituent = static_cast<Candidate *>(iterator->Next())))
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
