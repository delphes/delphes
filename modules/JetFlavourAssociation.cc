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
 *  MERCHANTABILITY || FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/** \class JetFlavourAssociation
 *
 *  Find origin of jet && evaluate jet flavor
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/JetFlavourAssociation.h"

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

class PartonClassifier: public ExRootClassifier
{
public:

  PartonClassifier() {}
  Int_t GetCategory(TObject *object);
  Double_t fEtaMax, fPTMin;
};

//------------------------------------------------------------------------------
// https://cmssdt.cern.ch/SDT/lxr/source/PhysicsTools/JetMCAlgos/plugins/PartonSelector.cc

Int_t PartonClassifier::GetCategory(TObject *object)
{
 // select parton in the parton list

  Candidate *parton = static_cast<Candidate *>(object);
  const TLorentzVector &momentum = parton->Momentum;
  Int_t pdgCode;

  // inside the eta && momentum range (be a little bit larger that the tracking coverage
  if(momentum.Pt() <= fPTMin || TMath::Abs(momentum.Eta()) > fEtaMax) return -1;

  pdgCode = TMath::Abs(parton->PID);

  if(parton->Status == -1) return -1;
  if(pdgCode != 21 && pdgCode > 5) return -1; // not a parton, skip
  if(parton->Status == 3 || parton->Status == 2) return 0; // if status 3 return

  return 0;
}

//------------------------------------------------------------------------------

class PartonClassifierLHEF: public ExRootClassifier
{
public:

  PartonClassifierLHEF() {}
  Int_t GetCategory(TObject *object);
  Double_t fEtaMax, fPTMin;
};

Int_t PartonClassifierLHEF::GetCategory(TObject *object)
{
  // select parton in the parton list

  Candidate *parton = static_cast<Candidate *>(object);
  const TLorentzVector &momentum = parton->Momentum;
  Int_t pdgCode;

  // inside the eta && momentum range (be a little bit larger that the tracking coverage
  if(momentum.Pt() <= fPTMin || TMath::Abs(momentum.Eta()) > fEtaMax) return -1;

  pdgCode = TMath::Abs(parton->PID);
  if(parton->Status == -1) return -1;
  if(pdgCode != 21 && pdgCode > 5) return -1; // not a parton, skip
  if(parton->Status != 1) return -1; // if status 3 return

  return 0;
}

//------------------------------------------------------------------------------

JetFlavourAssociation::JetFlavourAssociation() :
  fClassifier(0), fFilter(0),
  fItPartonInputArray(0), fItPartonInputArrayLHEF(0),
  fItJetInputArray(0), fItParticleInputArray(0)
{
  fClassifier    = new PartonClassifier;
  fClassifierLHEF = new PartonClassifierLHEF;
}

//------------------------------------------------------------------------------

JetFlavourAssociation::~JetFlavourAssociation()
{
  if(fClassifier) delete fClassifier;
  if(fClassifierLHEF) delete fClassifierLHEF;
}

//------------------------------------------------------------------------------

void JetFlavourAssociation::Init()
{
  ExRootConfParam param;

  fDeltaR = GetDouble("DeltaR", 0.5);

  fClassifier->fPTMin = GetDouble("PartonPTMin", 0.);
  fClassifier->fEtaMax = GetDouble("PartonEtaMax",2.5);

  fClassifierLHEF->fPTMin = GetDouble("PartonPTMin", 0.);
  fClassifierLHEF->fEtaMax = GetDouble("PartonEtaMax",2.5);

  // import input array(s)
  fPartonInputArray = ImportArray(GetString("PartonInputArray", "Delphes/partons"));
  fItPartonInputArray = fPartonInputArray->MakeIterator();

  fPartonInputArrayLHEF = ImportArray(GetString("PartonInputArrayLHEF", "Delphes/partonsLHEF"));
  fItPartonInputArrayLHEF = fPartonInputArrayLHEF->MakeIterator();

  fFilter = new ExRootFilter(fPartonInputArray);
  fFilterLHEF = new ExRootFilter(fPartonInputArrayLHEF);

  fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();

  fParticleInputArray = ImportArray(GetString("ParticleInputArray", "Delphes/allParticles"));
  fItParticleInputArray = fParticleInputArray->MakeIterator();

}

//------------------------------------------------------------------------------

void JetFlavourAssociation::Finish()
{
  if(fFilter) delete fFilter;
  if(fFilterLHEF) delete fFilterLHEF;

  if(fItJetInputArray) delete fItJetInputArray;
  if(fItParticleInputArray) delete fItParticleInputArray;
  if(fItPartonInputArray) delete fItPartonInputArray;
  if(fItPartonInputArrayLHEF) delete fItPartonInputArrayLHEF;
}

//------------------------------------------------------------------------------

void JetFlavourAssociation::Process(){

  Candidate *jet;
  TObjArray *partonArray;
  TObjArray *partonArrayLHEF;

  // select quark && gluons
  fFilter->Reset();
  partonArray = fFilter->GetSubArray(fClassifier, 0); // get the filtered parton array

  if(partonArray == 0) return;
  TIter itPartonArray(partonArray);

  fFilterLHEF->Reset();
  partonArrayLHEF = fFilterLHEF->GetSubArray(fClassifierLHEF, 0); // get the filtered parton array

  if(partonArrayLHEF == 0) return;
  TIter itPartonLHEFArray(partonArrayLHEF);

  // loop over all input jets
  fItJetInputArray->Reset();
  while((jet = static_cast<Candidate *>(fItJetInputArray->Next())))
  {
    // get standard flavor
    GetAlgoFlavour(jet, itPartonArray, itPartonLHEFArray);
    GetPhysicsFlavour(jet, itPartonArray, itPartonLHEFArray);
  }
}

//------------------------------------------------------------------------------
// Standard definition of jet flavor in
// https://cmssdt.cern.ch/SDT/lxr/source/PhysicsTools/JetMCAlgos/plugins/JetPartonMatcher.cc?v=CMSSW_7_3_0_pre1

void JetFlavourAssociation::GetAlgoFlavour(Candidate *jet, TIter &itPartonArray, TIter &itPartonLHEFArray)
{
  float maxPt = 0;
  float minDr = 1000;
  bool isGoodParton = true;
  int daughterCounter = 0;
  Candidate *parton, *partonLHEF;
  Candidate *tempParton = 0, *tempPartonNearest = 0, *tempPartonHighestPt = 0;
  int pdgCode, pdgCodeMax = -1;

  itPartonArray.Reset();
  while((parton = static_cast<Candidate *>(itPartonArray.Next())))
  {
    // default delphes method
    pdgCode = TMath::Abs(parton->PID);
    if(TMath::Abs(parton->PID) == 21) pdgCode = 0;
    if(jet->Momentum.DeltaR(parton->Momentum) <= fDeltaR)
    {
      if(pdgCodeMax < pdgCode) pdgCodeMax = pdgCode;
    }

    isGoodParton = true;
    itPartonLHEFArray.Reset();
    while((partonLHEF = static_cast<Candidate *>(itPartonLHEFArray.Next())))
    {
      if(parton->Momentum.DeltaR(partonLHEF->Momentum) < 0.001 &&
         parton->PID == partonLHEF->PID &&
         partonLHEF->Charge == parton->Charge)
      {
         isGoodParton = false;
         break;
      }

      if(!isGoodParton) continue;

      // check the daugheter
      daughterCounter = 0;
      if(parton->D1 != -1 || parton->D2 != -1)
      {
        // partons are only quarks || gluons
        int daughterFlavour1 = -1;
        int daughterFlavour2 = -1;
        if(parton->D1 != -1) daughterFlavour1 = TMath::Abs(static_cast<Candidate *>(fParticleInputArray->At(parton->D1))->PID);
        if(parton->D2 != -1) daughterFlavour2 = TMath::Abs(static_cast<Candidate *>(fParticleInputArray->At(parton->D2))->PID);
        if((daughterFlavour1 == 1 || daughterFlavour1 == 2 || daughterFlavour1 == 3 || daughterFlavour1 == 4 || daughterFlavour1 == 5 || daughterFlavour1 == 21)) daughterCounter++;
        if((daughterFlavour2 == 1 || daughterFlavour2 == 2 || daughterFlavour2 == 3 || daughterFlavour2 == 4 || daughterFlavour1 == 5 || daughterFlavour2 == 21)) daughterCounter++;
      }
      if(daughterCounter > 0) continue;
      if(jet->Momentum.DeltaR(parton->Momentum) <= fDeltaR)
      {
        if(jet->Momentum.DeltaR(parton->Momentum) < minDr)
        {
          minDr = jet->Momentum.DeltaR(parton->Momentum);
          tempPartonNearest = parton;
        }

        // if not yet found && pdgId is a c, take as c
        if(TMath::Abs(parton->PID) == 4) tempParton = parton;
        if(TMath::Abs(parton->PID) == 5) tempParton = parton;
        if(parton->Momentum.Pt() > maxPt)
        {
          maxPt = parton->Momentum.Pt();
          tempPartonHighestPt = parton;
        }
      }
    }
  }

  jet->FlavorHeaviest = tempParton ? TMath::Abs(tempParton->PID) : 0;
  jet->FlavorHighestPt = tempPartonHighestPt ? TMath::Abs(tempPartonHighestPt->PID) : 0;
  jet->FlavorNearest2 = tempPartonNearest ? TMath::Abs(tempPartonNearest->PID) : 0;
  if(!tempParton) tempParton = tempPartonHighestPt;
  jet->FlavorAlgo = tempParton ? TMath::Abs(tempParton->PID) : 0;

  if(pdgCodeMax == 0) pdgCodeMax = 21;
  if(pdgCodeMax == -1) pdgCodeMax = 0;

  jet->FlavorDefault = pdgCodeMax;
}

//------------------------------------------------------------------------------

void JetFlavourAssociation::GetPhysicsFlavour(Candidate *jet, TIter &itPartonArray, TIter &itPartonLHEFArray)
{
  float minDr = 1000;
  int partonCounter = 0;
  float biggerConeSize = 0.7;
  float dist;
  bool isGoodCandidate;
  int contaminatingFlavour = 0;
  int motherCounter = 0;
  Candidate *parton, *partonLHEF, *mother1, *mother2;
  Candidate *tempParton = 0, *tempPartonNearest = 0;
  vector<Candidate *> contaminations;
  vector<Candidate *>::iterator itContaminations;

  contaminations.clear();

  itPartonLHEFArray.Reset();
  while((partonLHEF = static_cast<Candidate *>(itPartonLHEFArray.Next())))
  {
    dist = jet->Momentum.DeltaR(partonLHEF->Momentum); // take the DR
    if(partonLHEF->Status == 1 && dist < minDr)
    {
      tempPartonNearest = partonLHEF;
      minDr = dist;
    }

    if(partonLHEF->Status == 1 && dist <= fDeltaR)
    {
      tempParton = partonLHEF;
      partonCounter++;
    }
  }

  itPartonArray.Reset();
  itPartonLHEFArray.Reset();
  while((parton = static_cast<Candidate *>(itPartonArray.Next())))
  {
    dist = jet->Momentum.DeltaR(parton->Momentum); // take the DR
    isGoodCandidate = true;
    while((partonLHEF = static_cast<Candidate *>(itPartonLHEFArray.Next())))
    {
      if(parton->Momentum.DeltaR(partonLHEF->Momentum) < 0.01 &&
         parton->PID == partonLHEF->PID &&
         partonLHEF->Charge == parton->Charge)
      {
        isGoodCandidate = false;
        break;
      }
    }

    if(!isGoodCandidate) continue;

    if(parton->D1 != -1 || parton->D2 != -1)
    {
      if((TMath::Abs(parton->PID) < 4 || TMath::Abs(parton->PID) == 21)) continue;
      if(dist < biggerConeSize) contaminations.push_back(parton);
    }
  }

  jet->FlavorNearest3 = tempPartonNearest ? TMath::Abs(tempPartonNearest->PID) : 0;

  if(partonCounter != 1)
  {
    jet->FlavorPhysics = 0;
  }
  else if(contaminations.size() == 0)
  {
    jet->FlavorPhysics = TMath::Abs(tempParton->PID);
  }
  else if(contaminations.size() > 0)
  {
    jet->FlavorPhysics = TMath::Abs(tempParton->PID);

    for(itContaminations = contaminations.begin(); itContaminations != contaminations.end(); ++itContaminations)
    {
      parton = *itContaminations;
      contaminatingFlavour = TMath::Abs(parton->PID);
      motherCounter = 0;
      if(parton->M1 != -1) motherCounter++;
      if(parton->M2 != -1) motherCounter++;

      if(parton->M1 != -1)
      {
        mother1 = static_cast<Candidate *>(fParticleInputArray->At(parton->M1));
        if(mother1 && motherCounter > 0 && mother1->Momentum.DeltaR(tempParton->Momentum) < 0.001) continue;
      }
      if(parton->M2 != -1)
      {
        mother2 = static_cast<Candidate *>(fParticleInputArray->At(parton->M2));
        if(mother2 && motherCounter > 0 && mother2->Momentum.DeltaR(tempParton->Momentum) < 0.001) continue;
      }
      // mother is the initialParton --> OK
      if(TMath::Abs(tempParton->PID) == 4)
      {
        // keep association --> the initialParton is a c --> the contaminated parton is a c
        if(contaminatingFlavour == 4) continue;
        jet->FlavorPhysics = 0; // all the other cases reject!
        break;
      }
    }
  }
}
