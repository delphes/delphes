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

/** \class TauTagging
 *
 *  Determines origin of jet,
 *  applies b-tagging efficiency (miss identification rate) formulas
 *  and sets b-tagging flags
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootSTLVectorFilter.h"

#include "modules/TauTagging.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

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

using namespace std;

//------------------------------------------------------------------------------

Int_t TauTaggingPartonClassifier::GetCategory(TObject *object)
{
  Candidate *tau = static_cast<Candidate *>(object);

  const TLorentzVector &momentum = tau->Momentum;
  Int_t pdgCode;

  pdgCode = TMath::Abs(tau->PID);
  if(pdgCode != 15) return -1;

  if(momentum.Pt() <= fPTMin || TMath::Abs(momentum.Eta()) > fEtaMax) return -1;

  if(tau->D1 < 0) return -1;

  if(tau->D2 < tau->D1) return -1;

  if(tau->D1 >= static_cast<int>(fParticleInputArray.size()) || tau->D2 >= static_cast<int>(fParticleInputArray.size()))
  {
    throw runtime_error("tau's daughter index is greater than the ParticleInputArray size");
  }

  for(int i = tau->D1; i <= tau->D2; ++i)
  {
    const auto &daughter1 = fParticleInputArray.at(i);
    pdgCode = TMath::Abs(daughter1.PID);
    //if(pdgCode == 11 || pdgCode == 13 || pdgCode == 15)
    //  return -1;
    if(pdgCode == 24)
    {
      if(daughter1.D1 < 0) return -1;
      for(int j = daughter1.D1; j <= daughter1.D2; ++j)
      {
        const auto &daughter2 = fParticleInputArray.at(j);
        pdgCode = TMath::Abs(daughter2.PID);
        if(pdgCode == 11 || pdgCode == 13) return -1;
      }
    }
  }
  return 0;
}

//------------------------------------------------------------------------------

TauTagging::TauTagging() :
  fClassifier(0), fFilter(0)
{
}

//------------------------------------------------------------------------------

TauTagging::~TauTagging()
{
}

//------------------------------------------------------------------------------

void TauTagging::Init()
{
  map<Int_t, DelphesFormula *>::iterator itEfficiencyMap;
  ExRootConfParam param;
  DelphesFormula *formula;
  Int_t i, size;

  fBitNumber = GetInt("BitNumber", 0);

  fDeltaR = GetDouble("DeltaR", 0.5);

  // read efficiency formulas
  param = GetParam("EfficiencyFormula");
  size = param.GetSize();

  fEfficiencyMap.clear();
  for(i = 0; i < size / 2; ++i)
  {
    formula = new DelphesFormula;
    formula->Compile(param[i * 2 + 1].GetString());

    fEfficiencyMap[param[i * 2].GetInt()] = formula;
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
  GetFactory()->EventModel()->Attach(GetString("ParticleInputArray", "Delphes/allParticles"), fParticleInputArray);
  GetFactory()->EventModel()->Attach(GetString("PartonInputArray", "Delphes/partons"), fPartonInputArray);
  GetFactory()->EventModel()->Attach(GetString("JetInputArray", "FastJetFinder/jets"), fJetInputArray); // I/O

  fClassifier = new TauTaggingPartonClassifier(*fParticleInputArray);
  fClassifier->fPTMin = GetDouble("TauPTMin", 1.0);
  fClassifier->fEtaMax = GetDouble("TauEtaMax", 2.5);

  fFilter = new ExRootSTLVectorFilter(*fPartonInputArray);
}

//------------------------------------------------------------------------------

void TauTagging::Finish()
{
  map<Int_t, DelphesFormula *>::iterator itEfficiencyMap;
  DelphesFormula *formula;

  if(fFilter) delete fFilter;
  if(fClassifier) delete fClassifier;

  for(itEfficiencyMap = fEfficiencyMap.begin(); itEfficiencyMap != fEfficiencyMap.end(); ++itEfficiencyMap)
  {
    formula = itEfficiencyMap->second;
    if(formula) delete formula;
  }
}

//------------------------------------------------------------------------------

void TauTagging::Process()
{
  TLorentzVector tauMomentum;
  Double_t pt, eta, phi, e, eff;
  map<Int_t, DelphesFormula *>::iterator itEfficiencyMap;
  DelphesFormula *formula;
  Int_t pdgCode, charge;

  // select taus
  fFilter->Reset();
  const auto tauArray = fFilter->GetSubArray(fClassifier, 0);

  // loop over all input jets
  for(auto &jet : *fJetInputArray)
  {
    const TLorentzVector &jetMomentum = jet.Momentum;
    pdgCode = 0;
    charge = gRandom->Uniform() > 0.5 ? 1 : -1;
    eta = jetMomentum.Eta();
    phi = jetMomentum.Phi();
    pt = jetMomentum.Pt();
    e = jetMomentum.E();

    // loop over all input taus
    if(!tauArray.empty())
    {
      for(const auto &tau : tauArray)
      {
        if(tau.D1 < 0) continue;

        if(tau.D1 >= static_cast<int>(fParticleInputArray->size()) || tau.D2 >= static_cast<int>(fParticleInputArray->size()))
        {
          throw runtime_error("tau's daughter index is greater than the ParticleInputArray size");
        }

        tauMomentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);

        for(int i = tau.D1; i <= tau.D2; ++i)
        {
          const auto &daughter = fParticleInputArray->at(i);
          if(TMath::Abs(daughter.PID) == 16) continue;
          tauMomentum += daughter.Momentum;
        }

        if(jetMomentum.DeltaR(tauMomentum) <= fDeltaR)
        {
          pdgCode = 15;
          charge = tau.Charge;
        }
      }
    }

    // fake electrons and muons

    if(pdgCode == 0)
    {

      Double_t drMin = fDeltaR;
      for(const auto &part : *fPartonInputArray)
      {
        if(TMath::Abs(part.PID) == 11 || TMath::Abs(part.PID) == 13)
        {
          tauMomentum = part.Momentum;
          if(tauMomentum.Pt() < fClassifier->fPTMin) continue;
          if(TMath::Abs(tauMomentum.Eta()) > fClassifier->fEtaMax) continue;

          Double_t dr = jetMomentum.DeltaR(tauMomentum);
          if(dr < drMin)
          {
            drMin = dr;
            pdgCode = TMath::Abs(part.PID);
            charge = part.Charge;
          }
        }
      }
    }

    // find an efficency formula
    itEfficiencyMap = fEfficiencyMap.find(pdgCode);
    if(itEfficiencyMap == fEfficiencyMap.end())
    {
      itEfficiencyMap = fEfficiencyMap.find(0);
    }
    formula = itEfficiencyMap->second;

    // apply an efficency formula
    eff = formula->Eval(pt, eta, phi, e);
    jet.TauFlavor = pdgCode;
    jet.TauTag |= (gRandom->Uniform() <= eff) << fBitNumber;
    jet.TauWeight = eff;

    // set tau charge
    jet.Charge = charge;
  }
}

//------------------------------------------------------------------------------
