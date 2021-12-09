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

/** \class LLPFilter
 *
 *  Removes particles with specific PDG codes
 *
 *  \author M. Selvaggi
 *
 */

#include "modules/LLPFilter.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

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

LLPFilter::LLPFilter() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

LLPFilter::~LLPFilter()
{
}

//------------------------------------------------------------------------------

void LLPFilter::Init()
{

  ExRootConfParam param;
  Size_t i, size;

  // PT threshold
  fPTMin = GetDouble("PTMin", 0.0);

  fInvert = GetBool("Invert", false);
  fDecayRegion = GetInt("DecayRegion", 0);
  fDaughterNumber = GetInt("DaughterNumber", 0);



  // no pileup
  fRequireNotPileup = GetBool("RequireNotPileup", false);

  fRequireStatus = GetBool("RequireStatus", false);
  fStatus = GetInt("Status", 1);

  fRequireCharge = GetBool("RequireCharge", false);
  fCharge = GetInt("Charge", 1);

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));
  fItInputArray = fInputArray->MakeIterator();

  fParticleInputArray =  ImportArray(GetString("InputArray", "Delphes/allParticles"));
  fItParticleInputArray = fParticleInputArray->MakeIterator();


  param = GetParam("PdgCode");
  size = param.GetSize();

  // read PdgCodes to be filtered out from the data card

  fPdgCodes.clear();
  for(i = 0; i < size; ++i)
  {
    fPdgCodes.push_back(param[i].GetInt());
  }

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "filteredParticles"));
}

//------------------------------------------------------------------------------

void LLPFilter::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void LLPFilter::Process()
{

  Candidate *candidate;
  Int_t pdgCode;
  Bool_t pass;
  Double_t pt, eta;
  Candidate *tempCandidate;

  Candidate *daughter;
  Int_t daughterPdg;


  // loop over particles to find LLP
  fItInputArray->Reset();
  int index = -1;
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    index++;

    //all distance units are in mm
    pdgCode = candidate->PID;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    const TLorentzVector &candidateProdPosition = candidate->Position;
    const TLorentzVector &candidateDecayPosition = candidate->DecayPosition;
    pt = candidateMomentum.Pt();
    eta = candidateMomentum.Eta();
    if(pt < fPTMin) continue;
    if (fDaughterNumber > 0)
    {
      if (candidate->D2-candidate->D1 != fDaughterNumber) continue;//require 3 daughters

    }
    if(find(fPdgCodes.begin(), fPdgCodes.end(), pdgCode) == fPdgCodes.end()) continue; //require pdgID is one of the LLP id
    if(fRequireStatus && (candidate->Status != fStatus)) continue;

    // loop over particles to find LLP daughters and assign EM and hadronic energy
    candidate->Eem = 0.0;
    candidate->Ehad = 0.0;

    fItParticleInputArray->Reset();

    while((daughter = static_cast<Candidate *>(fItParticleInputArray->Next())))
    {

      // check mother is LLP
      //check ID is 123456 or neutrinos or charged leptons
      daughterPdg = daughter->PID;
      if (daughter->Status != 1)continue;

      if (daughter->IsPU)continue;

      if (abs(daughterPdg)==12 || abs(daughterPdg)==14 || abs(daughterPdg)==16 || abs(daughterPdg)==13)continue; // ignore neutrinos and muons


      if (abs(daughterPdg) > 1000000) continue;//ignore BSM particles
      const TLorentzVector &daughterProdPosition = daughter->Position;
      const TLorentzVector &daughterMomentum = daughter->Momentum;

      const TLorentzVector distance = daughterProdPosition - candidateDecayPosition;

      if (sqrt(pow(distance.X(), 2) + pow(distance.Y(), 2) + pow(distance.Z(), 2))>100) continue; //if vertices are close, then matched, this 1 m only works for those decay in muon system, would be too large for tracker related
      tempCandidate = daughter;
      while(tempCandidate->M1 != -1 && tempCandidate->M1 != index)
      {
        tempCandidate  = static_cast<Candidate *>(fParticleInputArray->At(tempCandidate->M1));
        // if (tempCandidate->PID == pdgCode)break;

      }
      if (tempCandidate->M1 == -1) continue;

      if (abs(daughterPdg)==11 || abs(daughterPdg)==22 || abs(daughterPdg)==111)candidate->Eem += daughterMomentum.E();
      else candidate->Ehad += daughterMomentum.E();


    }



    // used detector geometry in Figure 4.1.1, page141 from CERN-LHCC-97-032: https://cds.cern.ch/record/343814?ln=en
    // decayRegion = 0: no cuts on decay region
    // decayRegion = 1: select LLP that decays in CSC volume
    // decayRegion = 2: select LLP that decays outside of calorimeters
    if (fDecayRegion == 1)
    {
      if (abs(eta) < 2
         && abs(candidateDecayPosition.Z())<11000 && abs(candidateDecayPosition.Z())>4000
         && sqrt(pow(candidateDecayPosition.X(),2)+pow(candidateDecayPosition.Y(),2)) < 6955)
      {
        fOutputArray->Add(candidate);
      }

    }
    else if(fDecayRegion == 2)
    {
      if (abs(candidateDecayPosition.Z()) > 5680 && sqrt(pow(candidateDecayPosition.X(),2)+pow(candidateDecayPosition.Y(),2)) > 3000)
      {
        fOutputArray->Add(candidate);
      }
    }
    else{
      fOutputArray->Add(candidate);
    }
  }//end of while loop
}
