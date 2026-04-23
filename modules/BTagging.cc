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

#include <TRandom3.h>

#include <map>

class BTagging: public DelphesModule
{
public:
  explicit BTagging(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fBitNumber(Steer<int>("BitNumber", 0))
  {
    for(const std::pair<int, std::string> &efficiencyFormula :
      Steer<std::vector<std::pair<int, std::string> > >("EfficiencyFormula", {}))
    {
      std::unique_ptr<DelphesFormula> formula = std::make_unique<DelphesFormula>();
      formula->Compile(efficiencyFormula.second);
      fEfficiencyMap[efficiencyFormula.first] = std::move(formula);
    }
    if(fEfficiencyMap.count(0) == 0)
    { // set default efficiency formula
      std::unique_ptr<DelphesFormula> formula = std::make_unique<DelphesFormula>();
      formula->Compile("0.0");
      fEfficiencyMap[0] = std::move(formula);
    }
  }

  void Init() override
  {
    fJetInputArray = ImportArray(Steer<std::string>("JetInputArray", "FastJetFinder/jets"));
  }
  void Process() override;

private:
  const Int_t fBitNumber;

#if !defined(__CINT__) && !defined(__CLING__)
  std::map<Int_t, std::unique_ptr<DelphesFormula> > fEfficiencyMap; //!
#endif

  CandidatesCollection fJetInputArray; //!
};

using namespace std;

//------------------------------------------------------------------------------

void BTagging::Process()
{
  Double_t pt, eta, phi, e;
  std::map<Int_t, std::unique_ptr<DelphesFormula> >::iterator itEfficiencyMap;

  // loop over all input jets
  for(Candidate *const &jet : *fJetInputArray)
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
      std::unique_ptr<DelphesFormula> &formula = itEfficiencyMap->second;

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
      std::unique_ptr<DelphesFormula> &formula = itEfficiencyMap->second;

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
      std::unique_ptr<DelphesFormula> &formula = itEfficiencyMap->second;

      // apply an efficiency formula
      jet->BTagPhys |= (gRandom->Uniform() <= formula->Eval(pt, eta, phi, e)) << fBitNumber;
    }
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("BTagging", BTagging);
