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


/** \class Weighter
 *
 *  Apply a weight depending on PDG code.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/Weighter.h"

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

bool Weighter::TIndexStruct::operator< (const Weighter::TIndexStruct &value) const
{
  Int_t i;

  for(i = 0; i < 4; ++i)
  {
    if(codes[i] != value.codes[i]) return codes[i] < value.codes[i];
  }

  return false;
}

//------------------------------------------------------------------------------

Weighter::Weighter() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

Weighter::~Weighter()
{
}

//------------------------------------------------------------------------------

void Weighter::Init()
{
  ExRootConfParam param, paramCodes;
  Int_t i, j, size, sizeCodes;
  Int_t code;
  TIndexStruct index;
  Double_t weight;

  fWeightSet.clear();


  // set default weight value
  fWeightMap.clear();
  memset(index.codes, 0, 4*sizeof(Int_t));
  fWeightMap[index] = 1.0;

  // read weights
  param = GetParam("Weight");
  size = param.GetSize();
  for(i = 0; i < size/2; ++i)
  {
    paramCodes = param[i*2];
    sizeCodes = paramCodes.GetSize();
    weight = param[i*2 + 1].GetDouble();

    if(sizeCodes < 1 || sizeCodes > 4)
    {
      throw runtime_error("only 1, 2, 3 or 4 PDG codes can be specified per weight");
    }

    memset(index.codes, 0, 4*sizeof(Int_t));

    for(j = 0; j < sizeCodes; ++j)
    {
      code = paramCodes[j].GetInt();
      index.codes[j] = code;
      fWeightSet.insert(code);
    }

    sort(index.codes, index.codes + 4);

    fWeightMap[index] = weight;
  }

  // import input array(s)

  fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array(s)

  fOutputArray = ExportArray(GetString("OutputArray", "weight"));
}

//------------------------------------------------------------------------------

void Weighter::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void Weighter::Process()
{
  Candidate *candidate;
  Int_t i;
  TIndexStruct index;
  Double_t weight;
  set<Int_t>::iterator itCodeSet;
  map<TIndexStruct, Double_t>::iterator itWeightMap;

  DelphesFactory *factory = GetFactory();

  // loop over all particles
  fCodeSet.clear();
  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    if(candidate->Status != 3) continue;

    if(fWeightSet.find(candidate->PID) == fWeightSet.end()) continue;

    fCodeSet.insert(candidate->PID);
  }

  // find default weight value
  memset(index.codes, 0, 4*sizeof(Int_t));
  itWeightMap = fWeightMap.find(index);
  weight = itWeightMap->second;

  if(fCodeSet.size() <= 4)
  {
    i = 0;
    for(itCodeSet = fCodeSet.begin(); itCodeSet != fCodeSet.end(); ++itCodeSet)
    {
      index.codes[i] = *itCodeSet;
      ++i;
    }

    sort(index.codes, index.codes + 4);

    itWeightMap = fWeightMap.find(index);
    if(itWeightMap != fWeightMap.end())
    {
      weight = itWeightMap->second;
    }
  }

  candidate = factory->NewCandidate();
  candidate->Momentum.SetPtEtaPhiE(weight, 0.0, 0.0, weight);
  fOutputArray->Add(candidate);
}

//------------------------------------------------------------------------------
