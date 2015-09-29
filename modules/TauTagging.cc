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

#include "modules/TauTagging.h"

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
TauTaggingPartonClassifier::TauTaggingPartonClassifier(const TObjArray *array) :
  fParticleInputArray(array)
{
}

//------------------------------------------------------------------------------

Int_t TauTaggingPartonClassifier::GetCategory(TObject *object)
{
  Candidate *tau = static_cast<Candidate *>(object);
  Candidate *daughter1 = 0;
  Candidate *daughter2 = 0;
  
  const TLorentzVector &momentum = tau->Momentum;
  Int_t pdgCode, i, j;

  pdgCode = TMath::Abs(tau->PID);
  if(pdgCode != 15) return -1;

  if(momentum.Pt() <= fPTMin || TMath::Abs(momentum.Eta()) > fEtaMax) return -1;

  if(tau->D1 < 0) return -1;

  if(tau->D2 < tau->D1) return -1;

  if(tau->D1 >= fParticleInputArray->GetEntriesFast() ||
     tau->D2 >= fParticleInputArray->GetEntriesFast())
  {
    throw runtime_error("tau's daughter index is greater than the ParticleInputArray size");
  }

  for(i = tau->D1; i <= tau->D2; ++i)
  {
    daughter1 = static_cast<Candidate *>(fParticleInputArray->At(i));
    pdgCode = TMath::Abs(daughter1->PID);
    if(pdgCode == 11 || pdgCode == 13 || pdgCode == 15) return -1;
    else if(pdgCode == 24)
    {
     if(daughter1->D1 < 0) return -1;
     for(j = daughter1->D1; j <= daughter1->D2; ++j)
     {
       daughter2 = static_cast<Candidate*>(fParticleInputArray->At(j));
       pdgCode = TMath::Abs(daughter2->PID);
       if(pdgCode == 11 || pdgCode == 13) return -1;
     }
	
    }
  }

  return 0;
}

//------------------------------------------------------------------------------

TauTagging::TauTagging() :
  fClassifier(0), fFilter(0),
  fItPartonInputArray(0), fItJetInputArray(0)
{
}

//------------------------------------------------------------------------------

TauTagging::~TauTagging()
{
}

//------------------------------------------------------------------------------

void TauTagging::Init()
{
  map< Int_t, DelphesFormula * >::iterator itEfficiencyMap;
  ExRootConfParam param;
  DelphesFormula *formula;
  Int_t i, size;

  fDeltaR = GetDouble("DeltaR", 0.5);

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

  fParticleInputArray = ImportArray(GetString("ParticleInputArray", "Delphes/allParticles"));

  fClassifier = new TauTaggingPartonClassifier(fParticleInputArray);
  fClassifier->fPTMin = GetDouble("TauPTMin", 1.0);
  fClassifier->fEtaMax = GetDouble("TauEtaMax", 2.5);

  fPartonInputArray = ImportArray(GetString("PartonInputArray", "Delphes/partons"));
  fItPartonInputArray = fPartonInputArray->MakeIterator();

  fFilter = new ExRootFilter(fPartonInputArray);

  fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();
}

//------------------------------------------------------------------------------

void TauTagging::Finish()
{
  map< Int_t, DelphesFormula * >::iterator itEfficiencyMap;
  DelphesFormula *formula;

  if(fFilter) delete fFilter;
  if(fClassifier) delete fClassifier;
  if(fItJetInputArray) delete fItJetInputArray;
  if(fItPartonInputArray) delete fItPartonInputArray;

  for(itEfficiencyMap = fEfficiencyMap.begin(); itEfficiencyMap != fEfficiencyMap.end(); ++itEfficiencyMap)
  {
    formula = itEfficiencyMap->second;
    if(formula) delete formula;
  }
}

//------------------------------------------------------------------------------

void TauTagging::Process()
{
  Candidate *jet, *tau, *daughter;
  TLorentzVector tauMomentum;
  Double_t pt, eta, phi;
  TObjArray *tauArray;
  map< Int_t, DelphesFormula * >::iterator itEfficiencyMap;
  DelphesFormula *formula;
  Int_t pdgCode, charge, i;

  // select taus
  fFilter->Reset();
  tauArray = fFilter->GetSubArray(fClassifier, 0);

  // loop over all input jets
  fItJetInputArray->Reset();
  while((jet = static_cast<Candidate *>(fItJetInputArray->Next())))
  {
    const TLorentzVector &jetMomentum = jet->Momentum;
    pdgCode = 0;
    charge = gRandom->Uniform() > 0.5 ? 1 : -1;
    eta = jetMomentum.Eta();
    phi = jetMomentum.Phi();
    pt = jetMomentum.Pt();

    // loop over all input taus
    if(tauArray){
      TIter itTauArray(tauArray);
      while((tau = static_cast<Candidate *>(itTauArray.Next())))
      {
        if(tau->D1 < 0) continue;

        if(tau->D1 >= fParticleInputArray->GetEntriesFast() ||
           tau->D2 >= fParticleInputArray->GetEntriesFast())
        {
          throw runtime_error("tau's daughter index is greater than the ParticleInputArray size");
        }

        tauMomentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);

        for(i = tau->D1; i <= tau->D2; ++i)
        {
          daughter = static_cast<Candidate *>(fParticleInputArray->At(i));
          if(TMath::Abs(daughter->PID) == 16) continue;
          tauMomentum += daughter->Momentum;
        }

        if(jetMomentum.DeltaR(tauMomentum) <= fDeltaR)
        {
          pdgCode = 15;
          charge = tau->Charge;
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
    jet->TauTag = gRandom->Uniform() <= formula->Eval(pt, eta);
    // set tau charge
    jet->Charge = charge;
  }
}

//------------------------------------------------------------------------------
