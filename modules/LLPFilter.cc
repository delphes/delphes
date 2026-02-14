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
 *  Filter LLPs with particular PDG ID/status and calculate the EM and hadronic energy of LLP based on decay particles
 *  The classification of EM and hadronic energy of LLP is based on instructions from the HEPData entry for the CMS paper searching
 *  for neutral LLPs in the CMS endcap muon detectors: https://www.hepdata.net/record/104408
 *  Muons and neutrinos are ignored. Photons, electrons, and pi0 are EM energy and everything else is hadronic energy.
 *
 *  \author Christina Wang
 *
 */

#include "modules/LLPFilter.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

void LLPFilter::Init()
{

  ExRootConfParam param;
  Size_t i, size;

  // PT threshold
  fPTMin = GetDouble("PTMin", 0.0);

  fInvert = GetBool("Invert", false);
  fDaughterNumber = GetInt("DaughterNumber", 0);
  fRequireDecayRegion = GetBool("RequireDecayRegion", 0);

  fDecayRegionRMax = GetDouble("DecayRegionRMax", 0.0); //mm
  fDecayRegionRMin = GetDouble("DecayRegionRMin", 0.0); //mm
  fDecayRegionZMax = GetDouble("DecayRegionZMax", 0.0); //mm
  fDecayRegionZMin = GetDouble("DecayRegionZMin", 0.0); //mm
  fDecayRegionEtaMax = GetDouble("DecayRegionEtaMax", 0.0); // requirement on abs(eta)
  fDecayRegionEtaMin = GetDouble("DecayRegionEtaMin", 0.0); //requirement on abs(eta)

  // no pileup
  fRequireNotPileup = GetBool("RequireNotPileup", false);

  fRequireStatus = GetBool("RequireStatus", false);
  fStatus = GetInt("Status", 1);

  fRequireCharge = GetBool("RequireCharge", false);
  fCharge = GetInt("Charge", 1);

  // import input arrays
  ImportArray(GetString("InputArray", "Delphes/allParticles"), fInputArray); // I/O
  ImportArray(GetString("InputArray", "Delphes/allParticles"), fParticleInputArray);

  param = GetParam("PdgCode");
  size = param.GetSize();

  // read PdgCodes to be filtered out from the data card

  fPdgCodes.clear();
  for(i = 0; i < size; ++i)
  {
    fPdgCodes.push_back(param[i].GetInt());
  }

  // create output array
  ExportArray(fOutputArray, GetString("OutputArray", "filteredParticles"));
}

//------------------------------------------------------------------------------

void LLPFilter::Finish()
{
}

//------------------------------------------------------------------------------

void LLPFilter::Process()
{
  Int_t pdgCode;
  Double_t pt, eta;
  const Candidate *tempCandidate{nullptr};

  Int_t daughterPdg;

  // loop over particles to find LLP
  int index = -1;
  for(auto &candidate : *fInputArray)
  {
    index++;

    //all distance units are in mm
    pdgCode = candidate.PID;
    const auto &candidateMomentum = candidate.Momentum;
    const auto &candidateDecayPosition = candidate.DecayPosition;
    pt = candidateMomentum.Pt();
    eta = candidateMomentum.Eta();
    if(pt < fPTMin) continue;
    if(fDaughterNumber > 0)
    {
      if(candidate.D2 - candidate.D1 != fDaughterNumber) continue; //require at least fDaughterNumber daughters
    }
    if(find(fPdgCodes.begin(), fPdgCodes.end(), pdgCode) == fPdgCodes.end()) continue; //require pdgID is one of the LLP id
    if(fRequireStatus && (candidate.Status != fStatus)) continue;

    // loop over particles to find LLP daughters and assign EM and hadronic energy
    candidate.Eem = 0.0;
    candidate.Ehad = 0.0;

    for(const auto &daughter : *fParticleInputArray)
    {
      daughterPdg = daughter.PID;
      if(daughter.Status != 1) continue;
      if(daughter.IsPU) continue;
      if(abs(daughterPdg) == 12 || abs(daughterPdg) == 14 || abs(daughterPdg) == 16 || abs(daughterPdg) == 13) continue; // ignore neutrinos and muons
      if(abs(daughterPdg) > 1000000) continue; //ignore BSM particles

      const auto &daughterMomentum = daughter.Momentum;

      // look for mother until find LLP or reach the top of the tree
      tempCandidate = &daughter;
      while(tempCandidate->M1 != -1 && tempCandidate->M1 != index)
      {
        tempCandidate = static_cast<Candidate *>(&fParticleInputArray->at(tempCandidate->M1));
      }
      if(tempCandidate->M1 == -1) continue;

      // assign LLP EM or hadronic energy, depending on the daughter ID
      if(abs(daughterPdg) == 11 || abs(daughterPdg) == 22 || abs(daughterPdg) == 111)
        candidate.Eem += daughterMomentum.E();
      else
        candidate.Ehad += daughterMomentum.E();
    }

    if(fRequireDecayRegion)
    {
      if(abs(eta) < fDecayRegionEtaMax && abs(eta) > fDecayRegionEtaMin
        && abs(candidateDecayPosition.Z()) < fDecayRegionZMax && abs(candidateDecayPosition.Z()) > fDecayRegionZMin
        && sqrt(pow(candidateDecayPosition.X(), 2) + pow(candidateDecayPosition.Y(), 2)) < fDecayRegionRMax
        && sqrt(pow(candidateDecayPosition.X(), 2) + pow(candidateDecayPosition.Y(), 2)) > fDecayRegionRMin)
      {
        fOutputArray->emplace_back(candidate);
      }
    }
    else
    {
      fOutputArray->emplace_back(candidate);
    }
  } //end of while loop
}
