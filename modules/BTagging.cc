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


/** \class BTagging
 *
 *  Determines origin of jet,
 *  applies b-tagging efficiency (miss identification rate) formulas
 *  and sets b-tagging flags
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/BTagging.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

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

BTagging::BTagging() :
  fItJetInputArray(0)
{
}

//------------------------------------------------------------------------------

BTagging::~BTagging()
{
}

//------------------------------------------------------------------------------

void BTagging::Init()
{
  map< Int_t, DelphesFormula * >::iterator itEfficiencyMap;
  ExRootConfParam param;
  DelphesFormula *formula;
  Int_t i, size;

  fBitNumber = GetInt("BitNumber", 0);

  // read efficiency formulas
  param = GetParam("EfficiencyFormula");
  size = param.GetSize();

  fEfficiencyMap.clear();
  for(i = 0; i < size/2; ++i)
  {
    formula = new DelphesFormula;
    formula->Compile(param[i*2 + 1].GetString());

    fEfficiencyMap[param[i*2].GetInt()] = formula;
  }

  // set default efficiency formula
  itEfficiencyMap = fEfficiencyMap.find(0);
  if(itEfficiencyMap == fEfficiencyMap.end())
  {
    formula = new DelphesFormula;
    formula->Compile("0.0");

    fEfficiencyMap[0] = formula;
  }

  // import input array(s)

  fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();
}

//------------------------------------------------------------------------------

void BTagging::Finish()
{
  map< Int_t, DelphesFormula * >::iterator itEfficiencyMap;
  DelphesFormula *formula;

  if(fItJetInputArray) delete fItJetInputArray;

  for(itEfficiencyMap = fEfficiencyMap.begin(); itEfficiencyMap != fEfficiencyMap.end(); ++itEfficiencyMap)
  {
    formula = itEfficiencyMap->second;
    if(formula) delete formula;
  }
}

//------------------------------------------------------------------------------

void BTagging::Process()
{
  Candidate *jet;
  Double_t pt, eta, phi, e;
  map< Int_t, DelphesFormula * >::iterator itEfficiencyMap;
  DelphesFormula *formula;

  // loop over all input jets
  fItJetInputArray->Reset();
  while((jet = static_cast<Candidate*>(fItJetInputArray->Next())))
  {
    const TLorentzVector &jetMomentum = jet->Momentum;
    eta = jetMomentum.Eta();
    phi = jetMomentum.Phi();
    pt = jetMomentum.Pt();
    e = jetMomentum.E();

    // find an efficiency formula
    itEfficiencyMap = fEfficiencyMap.find(jet->Flavor);
    if(itEfficiencyMap == fEfficiencyMap.end())
    {
      itEfficiencyMap = fEfficiencyMap.find(0);
    }
    formula = itEfficiencyMap->second;

    // apply an efficiency formula
    jet->BTag |= (gRandom->Uniform() <= formula->Eval(pt, eta, phi, e)) << fBitNumber;

    // find an efficiency formula for algo flavor definition
    itEfficiencyMap = fEfficiencyMap.find(jet->FlavorAlgo);
    if(itEfficiencyMap == fEfficiencyMap.end())
    {
      itEfficiencyMap = fEfficiencyMap.find(0);
    }
    formula = itEfficiencyMap->second;

    // apply an efficiency formula
    jet->BTagAlgo |= (gRandom->Uniform() <= formula->Eval(pt, eta, phi, e)) << fBitNumber;

    // find an efficiency formula for phys flavor definition
    itEfficiencyMap = fEfficiencyMap.find(jet->FlavorPhys);
    if(itEfficiencyMap == fEfficiencyMap.end())
    {
      itEfficiencyMap = fEfficiencyMap.find(0);
    }
    formula = itEfficiencyMap->second;

    // apply an efficiency formula
    jet->BTagPhys |= (gRandom->Uniform() <= formula->Eval(pt, eta, phi, e)) << fBitNumber;
  }
}

//------------------------------------------------------------------------------
