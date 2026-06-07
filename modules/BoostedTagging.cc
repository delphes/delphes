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

/** \class BoostedTagging
 *
 *  Determines the boosted origin of a large-R jet (W/Z/H/top),
 *  applies tagging efficiency (mis-identification rate) formulas
 *  and sets boosted-tagging flags
 *
 *  \author Sihyun Jeon - Boston University
 *
 */

#include "modules/BoostedTagging.h"

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

BoostedTaggingClassifier::BoostedTaggingClassifier(const TObjArray *array) :
  fParticleInputArray(array)
{
}

//------------------------------------------------------------------------------

// select last-copy W/Z/H/top resonances in the (pt, eta) acceptance
Int_t BoostedTaggingClassifier::GetCategory(TObject *object)
{
  Candidate *resonance = static_cast<Candidate *>(object);
  Int_t pdgCode = TMath::Abs(resonance->PID);

  if(pdgCode != 24 && pdgCode != 23 && pdgCode != 25 && pdgCode != 6) return -1;

  const TLorentzVector &momentum = resonance->Momentum;
  if(momentum.Pt() <= fPTMin || TMath::Abs(momentum.Eta()) > fEtaMax) return -1;

  // keep only the last copy (no daughter with the same PID)
  Int_t d1 = resonance->D1, d2 = resonance->D2, n = fParticleInputArray->GetEntriesFast();
  if(d1 >= 0 && d2 >= d1 && d1 < n)
  {
    if(d2 >= n) d2 = n - 1;
    for(Int_t i = d1; i <= d2; ++i)
      if(static_cast<Candidate *>(fParticleInputArray->At(i))->PID == resonance->PID) return -1;
  }

  return 0;
}

//------------------------------------------------------------------------------

BoostedTagging::BoostedTagging()
{
}

//------------------------------------------------------------------------------

BoostedTagging::~BoostedTagging()
{
}

//------------------------------------------------------------------------------

void BoostedTagging::Init()
{
  map<Int_t, DelphesFormula *>::iterator itEfficiencyMap;
  ExRootConfParam param;
  DelphesFormula *formula;
  Int_t i, size;

  fBitNumber = GetInt("BitNumber", 0);

  // matching cone: defaults to half the jet radius (R/2) when DeltaR is not set
  fJetRadius = GetDouble("JetRadius", 0.8);
  fDeltaR = GetDouble("DeltaR", 0.5 * fJetRadius);

  fJetPTMin = GetDouble("JetPTMin", 0.0);

  // optional soft-drop mass window on the jet (default: no cut)
  fSoftDropMassMin = GetDouble("SoftDropMassMin", 0.0);
  fSoftDropMassMax = GetDouble("SoftDropMassMax", 1.0e9);

  // read efficiency formulas (keyed by matched resonance |PDG|, 0 = default)
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
  fParticleInputArray = ImportArray(GetString("ParticleInputArray", "Delphes/allParticles"));

  fClassifier = new BoostedTaggingClassifier(fParticleInputArray);
  fClassifier->fPTMin = GetDouble("ResonancePTMin", 0.0);
  fClassifier->fEtaMax = GetDouble("ResonanceEtaMax", 5.0);

  fFilter = new ExRootFilter(fParticleInputArray);

  fJetInputArray = ImportArray(GetString("JetInputArray", "FatJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();
}

//------------------------------------------------------------------------------

void BoostedTagging::Finish()
{
  map<Int_t, DelphesFormula *>::iterator itEfficiencyMap;
  DelphesFormula *formula;

  if(fFilter) delete fFilter;
  if(fClassifier) delete fClassifier;
  if(fItJetInputArray) delete fItJetInputArray;

  for(itEfficiencyMap = fEfficiencyMap.begin(); itEfficiencyMap != fEfficiencyMap.end(); ++itEfficiencyMap)
  {
    formula = itEfficiencyMap->second;
    if(formula) delete formula;
  }
}

//------------------------------------------------------------------------------

void BoostedTagging::Process()
{
  Candidate *jet, *resonance;
  Double_t pt, eta, phi, e, dr, bestDR, eff, mass;
  Int_t pdgCode, origin;
  TObjArray *resonanceArray;
  map<Int_t, DelphesFormula *>::iterator itEfficiencyMap;
  DelphesFormula *formula;

  // select resonances
  fFilter->Reset();
  resonanceArray = fFilter->GetSubArray(fClassifier, 0);

  // loop over all input jets
  fItJetInputArray->Reset();
  while((jet = static_cast<Candidate *>(fItJetInputArray->Next())))
  {
    const TLorentzVector &jetMomentum = jet->Momentum;
    pt = jetMomentum.Pt();
    eta = jetMomentum.Eta();
    phi = jetMomentum.Phi();
    e = jetMomentum.E();

    if(pt < fJetPTMin) continue;

    // optional soft-drop mass window
    mass = jet->SoftDroppedJet.M();
    if(mass < fSoftDropMassMin || mass > fSoftDropMassMax) continue;

    // match the jet to the nearest configured resonance within DeltaR
    origin = 0;
    bestDR = fDeltaR;
    if(resonanceArray)
    {
      TIter itResonanceArray(resonanceArray);
      while((resonance = static_cast<Candidate *>(itResonanceArray.Next())))
      {
        pdgCode = TMath::Abs(resonance->PID);
        if(fEfficiencyMap.find(pdgCode) == fEfficiencyMap.end()) continue;

        dr = jetMomentum.DeltaR(resonance->Momentum);
        if(dr < bestDR)
        {
          bestDR = dr;
          origin = pdgCode;
        }
      }
    }

    // find an efficiency formula
    itEfficiencyMap = fEfficiencyMap.find(origin);
    if(itEfficiencyMap == fEfficiencyMap.end())
    {
      itEfficiencyMap = fEfficiencyMap.find(0);
    }
    formula = itEfficiencyMap->second;

    // apply an efficiency formula
    eff = formula->Eval(pt, eta, phi, e);
    jet->BoostedTag |= (gRandom->Uniform() <= eff) << fBitNumber;
  }
}

//------------------------------------------------------------------------------
