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

//------------------------------------------------------------------------------

/** \class LLPFilter
 *
 *  Filter LLPs with particular PDG ID/status and calculate the EM and hadronic energy of LLP based on decay particles
 *  The classification of EM and hadronic energy of LLP is based on instructions from the HEPData entry for the CMS paper searching
 *  for neutral LLPs in the CMS endcap muon detectors: https://www.hepdata.net/record/104408
 *
 *  \author Christina Wang
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>

#include <algorithm>
#include <vector>

class LLPFilter: public DelphesModule
{
public:
  explicit LLPFilter(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fPTMin(Steer<double>("PTMin", 0.0)), // PT threshold
    fInvert(Steer<bool>("Invert", false)),
    fDaughterNumber(Steer<int>("DaughterNumber", 0)),
    fRequireDecayRegion(Steer<bool>("RequireDecayRegion", 0)),
    fDecayRegionRMax(Steer<double>("DecayRegionRMax", 0.0)), //mm
    fDecayRegionRMin(Steer<double>("DecayRegionRMin", 0.0)), //mm
    fDecayRegionZMax(Steer<double>("DecayRegionZMax", 0.0)), //mm
    fDecayRegionZMin(Steer<double>("DecayRegionZMin", 0.0)), //mm
    fDecayRegionEtaMax(Steer<double>("DecayRegionEtaMax", 0.0)), // requirement on abs(eta)
    fDecayRegionEtaMin(Steer<double>("DecayRegionEtaMin", 0.0)), //requirement on abs(eta)
    fRequireNotPileup(Steer<bool>("RequireNotPileup", false)), // no pileup
    fRequireStatus(Steer<bool>("RequireStatus", false)),
    fStatus(Steer<int>("Status", 1)),
    fRequireCharge(Steer<bool>("RequireCharge", false)),
    fCharge(Steer<int>("Charge", 1)),
    fPdgCodes(Steer<std::vector<int> >("PdgCode")) // PdgCodes to be filtered out from the data card
  {
  }

  void Init() override
  {
    fInputArray = ImportArray(Steer<std::string>("InputArray", "Delphes/allParticles"));
    fParticleInputArray = ImportArray(Steer<std::string>("InputArray", "Delphes/allParticles"));
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "filteredParticles"));
  }
  void Process() override;

private:
  const double fPTMin; //!
  const bool fInvert; //!
  const int fDaughterNumber;
  const bool fRequireDecayRegion;
  const double fDecayRegionRMax;
  const double fDecayRegionRMin;
  const double fDecayRegionZMax;
  const double fDecayRegionZMin;
  const double fDecayRegionEtaMax;
  const double fDecayRegionEtaMin;
  const bool fRequireNotPileup; //!
  const bool fRequireStatus; //!
  const int fStatus; //!
  const bool fRequireCharge; //!
  const int fCharge; //!

  const std::vector<int> fPdgCodes;

  CandidatesCollection fInputArray; //!
  CandidatesCollection fParticleInputArray;

  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void LLPFilter::Process()
{
  fOutputArray->clear();

  // loop over particles to find LLP
  int index = -1;
  for(Candidate *const &candidate : *fInputArray)
  {
    index++;

    //all distance units are in mm
    const int pdgCode = candidate->PID;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    const TLorentzVector &candidateDecayPosition = candidate->DecayPosition;

    if(const double pt = candidateMomentum.Pt(); pt < fPTMin) continue;

    //require at least fDaughterNumber daughters
    if(fDaughterNumber > 0 && candidate->D2 - candidate->D1 != fDaughterNumber) continue;

    const double eta = candidateMomentum.Eta();
    if(std::find(fPdgCodes.begin(), fPdgCodes.end(), pdgCode) == fPdgCodes.end()) continue; //require pdgID is one of the LLP id
    if(fRequireStatus && (candidate->Status != fStatus)) continue;

    // loop over particles to find LLP daughters and assign EM and hadronic energy
    candidate->Eem = 0.0;
    candidate->Ehad = 0.0;

    for(Candidate *const &daughter : *fParticleInputArray)
    {
      const int daughterPdg = daughter->PID;
      if(daughter->Status != 1) continue;
      if(daughter->IsPU) continue;
      if(std::abs(daughterPdg) == 12 || std::abs(daughterPdg) == 14 || std::abs(daughterPdg) == 16 || std::abs(daughterPdg) == 13) continue; // ignore neutrinos and muons
      if(std::abs(daughterPdg) > 1000000) continue; //ignore BSM particles

      const TLorentzVector &daughterMomentum = daughter->Momentum;

      // look for mother until find LLP or reach the top of the tree
      Candidate *tempCandidate = daughter;
      while(tempCandidate->M1 != -1 && tempCandidate->M1 != index)
      {
        tempCandidate = static_cast<Candidate *>(fParticleInputArray->at(tempCandidate->M1));
      }
      if(tempCandidate->M1 == -1) continue;

      // assign LLP EM or hadronic energy, depending on the daughter ID
      if(std::abs(daughterPdg) == 11 || std::abs(daughterPdg) == 22 || std::abs(daughterPdg) == 111)
        candidate->Eem += daughterMomentum.E();
      else
        candidate->Ehad += daughterMomentum.E();
    }

    if(fRequireDecayRegion)
    {
      if(const double absEta = std::fabs(eta),
        candidateDecayPositionR = std::hypot(candidateDecayPosition.X(), candidateDecayPosition.Y());
        absEta < fDecayRegionEtaMax && absEta > fDecayRegionEtaMin
        && std::fabs(candidateDecayPosition.Z()) < fDecayRegionZMax && std::fabs(candidateDecayPosition.Z()) > fDecayRegionZMin
        && candidateDecayPositionR < fDecayRegionRMax && candidateDecayPositionR > fDecayRegionRMin)
        fOutputArray->emplace_back(candidate);
    }
    else
      fOutputArray->emplace_back(candidate);
  } //end of while loop
}

//------------------------------------------------------------------------------

REGISTER_MODULE("LLPFilter", LLPFilter);
