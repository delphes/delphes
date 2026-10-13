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

/** \class EventCounterFilter
 *
 *  Selects events by object multiplicity.
 *
 *  \author Sihyun Jeon (Boston University)
 *
 */

#include "modules/EventCounterFilter.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"

#include "classes/DelphesTreeWriter.h"

#include "ExRootAnalysis/ExRootConfReader.h"

#include "TObjArray.h"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace std;

//------------------------------------------------------------------------------

EventCounterFilter::EventCounterFilter() :
  fTotal(0), fPassed(0), fWriter(0)
{
}

//------------------------------------------------------------------------------

EventCounterFilter::~EventCounterFilter()
{
}

//------------------------------------------------------------------------------

void EventCounterFilter::Init()
{
  ExRootConfParam param = GetParam("Requirement");
  Long_t i, size;

  fInputArrays.clear();
  fOps.clear();
  fCounts.clear();
  fNames.clear();

  size = param.GetSize();
  if(size % 3 != 0)
  {
    throw runtime_error("EventCounterFilter: each Requirement needs 3 tokens: <InputArray> <op> <count>");
  }

  for(i = 0; i < size / 3; ++i)
  {
    string name = param[i * 3].GetString();
    string op = param[i * 3 + 1].GetString();
    Int_t count = param[i * 3 + 2].GetInt();

    fInputArrays.push_back(ImportArray(name.c_str()));
    fNames.push_back(name);
    fCounts.push_back(count);

    if(op == "eq" || op == "==")
      fOps.push_back(kEq);
    else if(op == "ne" || op == "!=")
      fOps.push_back(kNe);
    else if(op == "gt" || op == ">")
      fOps.push_back(kGt);
    else if(op == "ge" || op == ">=")
      fOps.push_back(kGe);
    else if(op == "lt" || op == "<")
      fOps.push_back(kLt);
    else if(op == "le" || op == "<=")
      fOps.push_back(kLe);
    else
    {
      stringstream message;
      message << "EventCounterFilter: unknown operator '" << op << "' (use eq ne gt ge lt le)";
      throw runtime_error(message.str());
    }
  }

  fTotal = 0;
  fPassed = 0;

  fWriter = dynamic_cast<DelphesTreeWriter *>(GetTreeWriter());
}

//------------------------------------------------------------------------------

void EventCounterFilter::Finish()
{
  cout << "** EventCounterFilter: kept " << fPassed << " / " << fTotal << " events" << endl;
}

//------------------------------------------------------------------------------

void EventCounterFilter::Process()
{
  Bool_t pass = kTRUE;
  size_t i;
  Int_t n;

  ++fTotal;

  for(i = 0; i < fInputArrays.size(); ++i)
  {
    n = fInputArrays[i]->GetEntriesFast();
    switch(fOps[i])
    {
    case kEq: pass = pass && (n == fCounts[i]); break;
    case kNe: pass = pass && (n != fCounts[i]); break;
    case kGt: pass = pass && (n > fCounts[i]); break;
    case kGe: pass = pass && (n >= fCounts[i]); break;
    case kLt: pass = pass && (n < fCounts[i]); break;
    case kLe: pass = pass && (n <= fCounts[i]); break;
    }
    if(!pass) break;
  }

  if(pass)
  {
    ++fPassed;
  }
  else if(fWriter)
  {
    fWriter->SetEventVeto();
  }
}

//------------------------------------------------------------------------------
