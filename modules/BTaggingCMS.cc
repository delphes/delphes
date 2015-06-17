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


/** \class BTaggingCMS
 *
 *  Applies b-tagging efficiency (miss identification rate) formulas
 *  and sets b-tagging flags
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/BTaggingCMS.h"

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

BTaggingCMS::BTaggingCMS() :
  fItJetInputArray(0)
{
}

//------------------------------------------------------------------------------

BTaggingCMS::~BTaggingCMS()
{
}

//------------------------------------------------------------------------------

void BTaggingCMS::Init()
{
  map< Int_t, DelphesFormula * >::iterator itEfficiencyMap;
  ExRootConfParam param;
  DelphesFormula *formula;
  Int_t i, size;

  // read efficiency formulas
  param = GetParam("EfficiencyFormulaLoose");
  size = param.GetSize();

  fEfficiencyMapLoose.clear();
  for(i = 0; i < size/2; ++i)
  {
    formula = new DelphesFormula;
    formula->Compile(param[i*2 + 1].GetString());
    fEfficiencyMapLoose[param[i*2].GetInt()] = formula;
  }

  // set default efficiency formula
  itEfficiencyMap = fEfficiencyMapLoose.find(0);
  if(itEfficiencyMap == fEfficiencyMapLoose.end())
  {
    formula = new DelphesFormula;
    formula->Compile("0.0");
    fEfficiencyMapLoose[0] = formula;
  }

  // read efficiency formulas
  param = GetParam("EfficiencyFormulaMedium");
  size = param.GetSize();

  fEfficiencyMapMedium.clear();
  for(i = 0; i < size/2; ++i)
  {
    formula = new DelphesFormula;
    formula->Compile(param[i*2 + 1].GetString());
    fEfficiencyMapMedium[param[i*2].GetInt()] = formula;
  }

  // set default efficiency formula
  itEfficiencyMap = fEfficiencyMapMedium.find(0);
  if(itEfficiencyMap == fEfficiencyMapMedium.end())
  {
    formula = new DelphesFormula;
    formula->Compile("0.0");
    fEfficiencyMapMedium[0] = formula;
  }

  // read efficiency formulas
  param = GetParam("EfficiencyFormulaTight");
  size = param.GetSize();

  fEfficiencyMapTight.clear();
  for(i = 0; i < size/2; ++i)
  {
    formula = new DelphesFormula;
    formula->Compile(param[i*2 + 1].GetString());
    fEfficiencyMapTight[param[i*2].GetInt()] = formula;
  }

  // set default efficiency formula
  itEfficiencyMap = fEfficiencyMapTight.find(0);
  if(itEfficiencyMap == fEfficiencyMapTight.end())
  {
    formula = new DelphesFormula;
    formula->Compile("0.0");
    fEfficiencyMapTight[0] = formula;
  }

  // import input array(s)
  fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();
}

//------------------------------------------------------------------------------

void BTaggingCMS::Finish()
{
  map< Int_t, DelphesFormula * >::iterator itEfficiencyMap;
  DelphesFormula *formula;

  if(fItJetInputArray) delete fItJetInputArray;

  for(itEfficiencyMap = fEfficiencyMapLoose.begin(); itEfficiencyMap != fEfficiencyMapLoose.end(); ++itEfficiencyMap)
  {
    formula = itEfficiencyMap->second;
    if(formula) delete formula;
  }

  for(itEfficiencyMap = fEfficiencyMapMedium.begin(); itEfficiencyMap != fEfficiencyMapMedium.end(); ++itEfficiencyMap)
  {
    formula = itEfficiencyMap->second;
    if(formula) delete formula;
  }

  for(itEfficiencyMap = fEfficiencyMapTight.begin(); itEfficiencyMap != fEfficiencyMapTight.end(); ++itEfficiencyMap)
  {
    formula = itEfficiencyMap->second;
    if(formula) delete formula;
  }
}

//------------------------------------------------------------------------------

void BTaggingCMS::Process()
{
  Candidate *jet;
  map< Int_t, DelphesFormula * >::iterator itEfficiencyMap;
  DelphesFormula *formula;
  float randomNumber;

  // loop over all input jets
  fItJetInputArray->Reset();

  while((jet = static_cast<Candidate*>(fItJetInputArray->Next())))
  {
    const TLorentzVector &jetMomentum = jet->Momentum;

    randomNumber = gRandom->Uniform();

    // --------------------------
    // Loose b-tagging
    // --------------------------

    // find heaviest flavor and b-tag
    if(fEfficiencyMapLoose.size() != 0)
    {
      itEfficiencyMap = fEfficiencyMapLoose.find(jet->FlavorHeaviest);
      if(itEfficiencyMap == fEfficiencyMapLoose.end())
      {
        itEfficiencyMap = fEfficiencyMapLoose.find(0);
      }
      formula = itEfficiencyMap->second;

      if(randomNumber <= formula->Eval(jetMomentum.Pt(), jetMomentum.Eta())) jet->BTagHeaviest = 1;

      // find hisghestpt flavor and b-tag
      itEfficiencyMap = fEfficiencyMapLoose.find(jet->FlavorHighestPt);
      if(itEfficiencyMap == fEfficiencyMapLoose.end())
      {
        itEfficiencyMap = fEfficiencyMapLoose.find(0);
      }
      formula = itEfficiencyMap->second;

      if(randomNumber <= formula->Eval(jetMomentum.Pt(), jetMomentum.Eta())) jet->BTagHighestPt = 1;

      // find nearest2 flavor and b-tag
      itEfficiencyMap = fEfficiencyMapLoose.find(jet->FlavorNearest2);
      if(itEfficiencyMap == fEfficiencyMapLoose.end())
      {
        itEfficiencyMap = fEfficiencyMapLoose.find(0);
      }
      formula = itEfficiencyMap->second;

      if(randomNumber <= formula->Eval(jetMomentum.Pt(), jetMomentum.Eta())) jet->BTagNearest2 = 1;

      // find nearest3 flavor and b-tag
      itEfficiencyMap = fEfficiencyMapLoose.find(jet->FlavorNearest3);
      if(itEfficiencyMap == fEfficiencyMapLoose.end())
      {
        itEfficiencyMap = fEfficiencyMapLoose.find(0);
      }
      formula = itEfficiencyMap->second;

      if(randomNumber <= formula->Eval(jetMomentum.Pt(), jetMomentum.Eta())) jet->BTagNearest3 = 1;

      // find nearest3 flavor and b-tag
      itEfficiencyMap = fEfficiencyMapLoose.find(jet->FlavorAlgo);
      if(itEfficiencyMap == fEfficiencyMapLoose.end())
      {
        itEfficiencyMap = fEfficiencyMapLoose.find(0);
      }
      formula = itEfficiencyMap->second;

      if(randomNumber <= formula->Eval(jetMomentum.Pt(), jetMomentum.Eta())) jet->BTagAlgo = 1;

      // find nearest3 flavor and b-tag
      itEfficiencyMap = fEfficiencyMapLoose.find(jet->FlavorPhysics);
      if(itEfficiencyMap == fEfficiencyMapLoose.end())
      {
        itEfficiencyMap = fEfficiencyMapLoose.find(0);
      }
      formula = itEfficiencyMap->second;

      if(randomNumber <= formula->Eval(jetMomentum.Pt(), jetMomentum.Eta())) jet->BTagPhysics = 1;

      // default
      itEfficiencyMap = fEfficiencyMapLoose.find(jet->FlavorDefault);
      if(itEfficiencyMap == fEfficiencyMapLoose.end())
      {
        itEfficiencyMap = fEfficiencyMapLoose.find(0);
      }
      formula = itEfficiencyMap->second;

      if(randomNumber <= formula->Eval(jetMomentum.Pt(), jetMomentum.Eta())) jet->BTagDefault = 1;
    }

    // --------------------------
    // Medium b-tagging
    // --------------------------

    // find heaviest flavor and b-tag
    if(fEfficiencyMapMedium.size() != 0)
    {
      itEfficiencyMap = fEfficiencyMapMedium.find(jet->FlavorHeaviest);
      if(itEfficiencyMap == fEfficiencyMapMedium.end())
      {
        itEfficiencyMap = fEfficiencyMapMedium.find(0);
      }
      formula = itEfficiencyMap->second;

      if(randomNumber <= formula->Eval(jetMomentum.Pt(), jetMomentum.Eta())) jet->BTagHeaviest = 2;

      // find hisghestpt flavor and b-tag
      itEfficiencyMap = fEfficiencyMapMedium.find(jet->FlavorHighestPt);
      if(itEfficiencyMap == fEfficiencyMapMedium.end())
      {
        itEfficiencyMap = fEfficiencyMapMedium.find(0);
      }
      formula = itEfficiencyMap->second;

      if(randomNumber <= formula->Eval(jetMomentum.Pt(), jetMomentum.Eta())) jet->BTagHighestPt = 2;

      // find nearest2 flavor and b-tag
      itEfficiencyMap = fEfficiencyMapMedium.find(jet->FlavorNearest2);
      if(itEfficiencyMap == fEfficiencyMapMedium.end())
      {
        itEfficiencyMap = fEfficiencyMapMedium.find(0);
      }
      formula = itEfficiencyMap->second;

      if(randomNumber <= formula->Eval(jetMomentum.Pt(), jetMomentum.Eta())) jet->BTagNearest2 = 2;

      // find nearest3 flavor and b-tag
      itEfficiencyMap = fEfficiencyMapMedium.find(jet->FlavorNearest3);
      if(itEfficiencyMap == fEfficiencyMapMedium.end())
      {
        itEfficiencyMap = fEfficiencyMapMedium.find(0);
      }
      formula = itEfficiencyMap->second;

      if(randomNumber <= formula->Eval(jetMomentum.Pt(), jetMomentum.Eta())) jet->BTagNearest3 = 2;

      // find nearest3 flavor and b-tag
      itEfficiencyMap = fEfficiencyMapMedium.find(jet->FlavorAlgo);
      if(itEfficiencyMap == fEfficiencyMapMedium.end())
      {
        itEfficiencyMap = fEfficiencyMapMedium.find(0);
      }
      formula = itEfficiencyMap->second;

      if(randomNumber <= formula->Eval(jetMomentum.Pt(), jetMomentum.Eta())) jet->BTagAlgo = 2;

      // find nearest3 flavor and b-tag
      itEfficiencyMap = fEfficiencyMapMedium.find(jet->FlavorPhysics);
      if(itEfficiencyMap == fEfficiencyMapMedium.end())
      {
        itEfficiencyMap = fEfficiencyMapMedium.find(0);
      }
      formula = itEfficiencyMap->second;

      if(randomNumber <= formula->Eval(jetMomentum.Pt(), jetMomentum.Eta())) jet->BTagPhysics = 2;

      // default
      itEfficiencyMap = fEfficiencyMapMedium.find(jet->FlavorDefault);
      if(itEfficiencyMap == fEfficiencyMapMedium.end())
      {
        itEfficiencyMap = fEfficiencyMapMedium.find(0);
      }
      formula = itEfficiencyMap->second;

      if(randomNumber <= formula->Eval(jetMomentum.Pt(), jetMomentum.Eta())) jet->BTagDefault = 2;
    }

    // --------------------------
    // Tight b-tagging
    // --------------------------

    // find heaviest flavor and b-tag
    if(fEfficiencyMapTight.size() != 0)
    {
      itEfficiencyMap = fEfficiencyMapTight.find(jet->FlavorHeaviest);
      if(itEfficiencyMap == fEfficiencyMapTight.end())
      {
        itEfficiencyMap = fEfficiencyMapTight.find(0);
      }
      formula = itEfficiencyMap->second;

      if(randomNumber <= formula->Eval(jetMomentum.Pt(), jetMomentum.Eta())) jet->BTagHeaviest = 3;

      // find hisghestpt flavor and b-tag
      itEfficiencyMap = fEfficiencyMapTight.find(jet->FlavorHighestPt);
      if(itEfficiencyMap == fEfficiencyMapTight.end())
      {
        itEfficiencyMap = fEfficiencyMapTight.find(0);
      }
      formula = itEfficiencyMap->second;

      if(randomNumber <= formula->Eval(jetMomentum.Pt(), jetMomentum.Eta())) jet->BTagHighestPt = 3;

      // find nearest2 flavor and b-tag
      itEfficiencyMap = fEfficiencyMapTight.find(jet->FlavorNearest2);
      if(itEfficiencyMap == fEfficiencyMapTight.end())
      {
        itEfficiencyMap = fEfficiencyMapTight.find(0);
      }
      formula = itEfficiencyMap->second;

      if(randomNumber <= formula->Eval(jetMomentum.Pt(), jetMomentum.Eta())) jet->BTagNearest2 = 3;

      // find nearest3 flavor and b-tag
      itEfficiencyMap = fEfficiencyMapTight.find(jet->FlavorNearest3);
      if(itEfficiencyMap == fEfficiencyMapTight.end())
      {
        itEfficiencyMap = fEfficiencyMapTight.find(0);
      }
      formula = itEfficiencyMap->second;

      if(randomNumber <= formula->Eval(jetMomentum.Pt(), jetMomentum.Eta())) jet->BTagNearest3 = 3;

      // find nearest3 flavor and b-tag
      itEfficiencyMap = fEfficiencyMapTight.find(jet->FlavorAlgo);
      if(itEfficiencyMap == fEfficiencyMapTight.end())
      {
        itEfficiencyMap = fEfficiencyMapTight.find(0);
      }
      formula = itEfficiencyMap->second;

      if(randomNumber <= formula->Eval(jetMomentum.Pt(), jetMomentum.Eta())) jet->BTagAlgo = 3;

      // find nearest3 flavor and b-tag
      itEfficiencyMap = fEfficiencyMapTight.find(jet->FlavorPhysics);
      if(itEfficiencyMap == fEfficiencyMapTight.end())
      {
        itEfficiencyMap = fEfficiencyMapTight.find(0);
      }
      formula = itEfficiencyMap->second;

      if(randomNumber <= formula->Eval(jetMomentum.Pt(), jetMomentum.Eta())) jet->BTagPhysics = 3;

      // default
      itEfficiencyMap = fEfficiencyMapTight.find(jet->FlavorDefault);
      if(itEfficiencyMap == fEfficiencyMapTight.end())
      {
        itEfficiencyMap = fEfficiencyMapTight.find(0);
      }
      formula = itEfficiencyMap->second;

      if(randomNumber <= formula->Eval(jetMomentum.Pt(), jetMomentum.Eta())) jet->BTagDefault = 3;
    }
  }
}
