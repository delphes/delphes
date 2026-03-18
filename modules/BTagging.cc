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

#include "classes/DelphesClasses.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"
#include "classes/DelphesModuleFactory.h"

#include <TRandom3.h>

#include <map>

class BTagging: public DelphesModule
{
public:
  BTagging() = default;

  void Init() override;
  void Process() override;

private:
  Int_t fBitNumber;

#if !defined(__CINT__) && !defined(__CLING__)
  std::map<Int_t, std::unique_ptr<DelphesFormula> > fEfficiencyMap; //!
#endif

  CandidatesCollection fJetInputArray; //!
};

using namespace std;

//------------------------------------------------------------------------------

void BTagging::Init()
{
  std::map<Int_t, std::unique_ptr<DelphesFormula> >::iterator itEfficiencyMap;
  ExRootConfParam param;
  Int_t i, size;

  fBitNumber = GetInt("BitNumber", 0);

  // read efficiency formulas
  param = GetParam("EfficiencyFormula");
  size = param.GetSize();

  fEfficiencyMap.clear();
  for(i = 0; i < size / 2; ++i)
  {
    auto formula = std::make_unique<DelphesFormula>();
    formula->Compile(param[i * 2 + 1].GetString());

    fEfficiencyMap[param[i * 2].GetInt()] = std::move(formula);
  }

  // set default efficiency formula
  itEfficiencyMap = fEfficiencyMap.find(0);
  if(itEfficiencyMap == fEfficiencyMap.end())
  {
    auto formula = std::make_unique<DelphesFormula>();
    formula->Compile("0.0");

    fEfficiencyMap[0] = std::move(formula);
  }

  // import input array(s)

  fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
}

//------------------------------------------------------------------------------

void BTagging::Process()
{
  Double_t pt, eta, phi, e;
  std::map<Int_t, std::unique_ptr<DelphesFormula> >::iterator itEfficiencyMap;

  // loop over all input jets
  for(const auto &jet : *fJetInputArray)
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
    {
      auto &formula = itEfficiencyMap->second;

      // apply an efficiency formula
      jet->BTag |= (gRandom->Uniform() <= formula->Eval(pt, eta, phi, e)) << fBitNumber;
    }

    // find an efficiency formula for algo flavor definition
    itEfficiencyMap = fEfficiencyMap.find(jet->FlavorAlgo);
    if(itEfficiencyMap == fEfficiencyMap.end())
    {
      itEfficiencyMap = fEfficiencyMap.find(0);
    }
    {
      auto &formula = itEfficiencyMap->second;

      // apply an efficiency formula
      jet->BTagAlgo |= (gRandom->Uniform() <= formula->Eval(pt, eta, phi, e)) << fBitNumber;
    }

    // find an efficiency formula for phys flavor definition
    itEfficiencyMap = fEfficiencyMap.find(jet->FlavorPhys);
    if(itEfficiencyMap == fEfficiencyMap.end())
    {
      itEfficiencyMap = fEfficiencyMap.find(0);
    }
    {
      auto &formula = itEfficiencyMap->second;

      // apply an efficiency formula
      jet->BTagPhys |= (gRandom->Uniform() <= formula->Eval(pt, eta, phi, e)) << fBitNumber;
    }
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("BTagging", BTagging);
