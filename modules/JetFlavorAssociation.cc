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

/** \class JetFlavorAssociation
 *
 *  Find origin of jet and evaluate jet flavor
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFilter.h"
#include "classes/DelphesModule.h"

#include "ExRootAnalysis/ExRootClassifier.h"

#include <TLorentzVector.h>

using namespace std;

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// https://cmssdt.cern.ch/SDT/lxr/source/PhysicsTools/JetMCAlgos/plugins/PartonSelector.cc
class PartonClassifier: public ExRootClassifier
{
public:
  PartonClassifier() {}
  int GetCategory(TObject *object)
  {
    // select parton in the parton list

    Candidate *parton = static_cast<Candidate *>(object);
    const TLorentzVector &momentum = parton->Momentum;
    int pdgCode;

    // inside the eta && momentum range (be a little bit larger that the tracking coverage
    if(momentum.Pt() <= fPTMin || std::fabs(momentum.Eta()) > fEtaMax) return -1;

    pdgCode = std::abs(parton->PID);

    if(parton->Status == -1) return -1;
    if(pdgCode != 21 && pdgCode > 5) return -1; // not a parton, skip
    if(parton->Status == 3 || parton->Status == 2) return 0; // if status 3 return

    return 0;
  }
  double fEtaMax, fPTMin;
};

//------------------------------------------------------------------------------

class ParticleLHEFClassifier: public ExRootClassifier
{
public:
  ParticleLHEFClassifier() {}
  int GetCategory(TObject *object)
  {
    // select parton in the parton list

    Candidate *particleLHEF = static_cast<Candidate *>(object);
    const TLorentzVector &momentum = particleLHEF->Momentum;
    int pdgCode;

    // inside the eta && momentum range (be a little bit larger that the tracking coverage
    if(momentum.Pt() <= fPTMin || std::fabs(momentum.Eta()) > fEtaMax) return -1;

    pdgCode = std::abs(particleLHEF->PID);
    if(particleLHEF->Status == -1) return -1;
    if(pdgCode != 21 && pdgCode > 5) return -1; // not a parton, skip
    if(particleLHEF->Status != 1) return -1; // if status 3 return

    return 0;
  }
  double fEtaMax, fPTMin;
};

//------------------------------------------------------------------------------

class JetFlavorAssociation: public DelphesModule
{
public:
  explicit JetFlavorAssociation(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fDeltaR(Steer<double>("DeltaR", 0.5)),
    fPartonFilter(std::make_unique<DelphesFilter>(fPartonInputArray)),
    fPartonClassifier(std::make_unique<PartonClassifier>()),
    fParticleLHEFClassifier(std::make_unique<ParticleLHEFClassifier>())
  {
    fPartonClassifier->fPTMin = Steer<double>("PartonPTMin", 0.0);
    fPartonClassifier->fEtaMax = Steer<double>("PartonEtaMax", 2.5);

    fParticleLHEFClassifier->fPTMin = Steer<double>("PartonPTMin", 0.0);
    fParticleLHEFClassifier->fEtaMax = Steer<double>("PartonEtaMax", 2.5);
  }

  void Init() override
  {
    fPartonInputArray = ImportArray(Steer<std::string>("PartonInputArray", "Delphes/partons"));
    fParticleInputArray = ImportArray(Steer<std::string>("ParticleInputArray", "Delphes/allParticles"));
    fJetInputArray = ImportArray(Steer<std::string>("JetInputArray", "FastJetFinder/jets"));
    try
    {
      fParticleLHEFInputArray = ImportArray(Steer<std::string>("ParticleLHEFInputArray", "Delphes/allParticlesLHEF"));
      fParticleLHEFFilter = std::make_unique<DelphesFilter>(fParticleLHEFInputArray);
    }
    catch(runtime_error &)
    {
    }
  }
  void Process() override;

  void GetAlgoFlavor(Candidate *jet, const std::vector<Candidate *> &partonArray, const std::vector<Candidate *> &partonLHEFArray);
  void GetPhysicsFlavor(Candidate *jet, const std::vector<Candidate *> &partonArray, const std::vector<Candidate *> &partonLHEFArray);

private:
  const double fDeltaR;

  const std::unique_ptr<DelphesFilter> fPartonFilter;
  const std::unique_ptr<PartonClassifier> fPartonClassifier; //!
  const std::unique_ptr<ParticleLHEFClassifier> fParticleLHEFClassifier; //!

  std::unique_ptr<DelphesFilter> fParticleLHEFFilter;

  CandidatesCollection fPartonInputArray; //!
  CandidatesCollection fParticleInputArray; //!
  CandidatesCollection fJetInputArray; //!
  CandidatesCollection fParticleLHEFInputArray; //!
};

//------------------------------------------------------------------------------

void JetFlavorAssociation::Process()
{

  // select quark and gluons
  fPartonFilter->Reset();
  const std::vector<Candidate *> partonArray = fPartonFilter->GetSubArray(fPartonClassifier.get(), 0); // get the filtered parton array
  if(partonArray.empty()) return;

  std::vector<Candidate *> partonLHEFArray;
  if(fParticleLHEFInputArray)
  {
    fParticleLHEFFilter->Reset();
    partonLHEFArray = fParticleLHEFFilter->GetSubArray(fParticleLHEFClassifier.get(), 0); // get the filtered parton array
  }
  // loop over all input jets
  for(Candidate *const &jet : *fJetInputArray)
  {
    // get standard flavor
    GetAlgoFlavor(jet, partonArray, partonLHEFArray);
    if(fParticleLHEFInputArray) GetPhysicsFlavor(jet, partonArray, partonLHEFArray);
  }
}

//------------------------------------------------------------------------------
// Standard definition of jet flavor in
// https://cmssdt.cern.ch/SDT/lxr/source/PhysicsTools/JetMCAlgos/plugins/JetPartonMatcher.cc?v=CMSSW_7_3_0_pre1

void JetFlavorAssociation::GetAlgoFlavor(Candidate *jet, const std::vector<Candidate *> &partonArray, const std::vector<Candidate *> &partonLHEFArray)
{
  float maxPt = 0;
  int daughterCounter = 0;
  Candidate *tempParton = 0, *tempPartonHighestPt = 0;
  int pdgCode, pdgCodeMax = -1;

  for(Candidate *const &parton : partonArray)
  {
    // default delphes method
    pdgCode = std::abs(parton->PID);
    if(std::abs(parton->PID) == 21) pdgCode = 0;
    if(jet->Momentum.DeltaR(parton->Momentum) <= fDeltaR)
    {
      if(pdgCodeMax < pdgCode) pdgCodeMax = pdgCode;
    }

    for(Candidate *const &partonLHEF : partonLHEFArray)
    {
      if(parton->Momentum.DeltaR(partonLHEF->Momentum) < 0.001 && parton->PID == partonLHEF->PID && partonLHEF->Charge == parton->Charge)
      {
        break;
      }

      // check the daughter
      daughterCounter = 0;
      if(parton->D1 != -1 || parton->D2 != -1)
      {
        // partons are only quarks || gluons
        int daughterFlavor1 = -1;
        int daughterFlavor2 = -1;
        if(parton->D1 != -1) daughterFlavor1 = std::abs(static_cast<Candidate *>(fParticleInputArray->at(parton->D1))->PID);
        if(parton->D2 != -1) daughterFlavor2 = std::abs(static_cast<Candidate *>(fParticleInputArray->at(parton->D2))->PID);
        if((daughterFlavor1 == 1 || daughterFlavor1 == 2 || daughterFlavor1 == 3 || daughterFlavor1 == 4 || daughterFlavor1 == 5 || daughterFlavor1 == 21)) daughterCounter++;
        if((daughterFlavor2 == 1 || daughterFlavor2 == 2 || daughterFlavor2 == 3 || daughterFlavor2 == 4 || daughterFlavor2 == 5 || daughterFlavor2 == 21)) daughterCounter++;
      }
      if(daughterCounter > 0) continue;
      if(jet->Momentum.DeltaR(parton->Momentum) <= fDeltaR)
      {
        // if not yet found && pdgId is a c, take as c
        if(std::abs(parton->PID) == 4) tempParton = parton;
        if(std::abs(parton->PID) == 5) tempParton = parton;
        if(parton->Momentum.Pt() > maxPt)
        {
          maxPt = parton->Momentum.Pt();
          tempPartonHighestPt = parton;
        }
      }
    }
  }

  if(!tempParton) tempParton = tempPartonHighestPt;
  jet->FlavorAlgo = tempParton ? std::abs(tempParton->PID) : 0;

  if(pdgCodeMax == 0) pdgCodeMax = 21;
  if(pdgCodeMax == -1) pdgCodeMax = 0;

  jet->Flavor = pdgCodeMax;
}

//------------------------------------------------------------------------------

void JetFlavorAssociation::GetPhysicsFlavor(Candidate *jet, const std::vector<Candidate *> &partonArray, const std::vector<Candidate *> &partonLHEFArray)
{
  int partonCounter = 0;
  float biggerConeSize = 0.7;
  float dist;
  bool isGoodCandidate;
  int contaminatingFlavor = 0;
  int motherCounter = 0;
  Candidate *tempParton = nullptr;
  vector<Candidate *> contaminations;

  contaminations.clear();

  for(Candidate *const &partonLHEF : partonLHEFArray)
  {
    dist = jet->Momentum.DeltaR(partonLHEF->Momentum); // take the DR

    if(partonLHEF->Status == 1 && dist <= fDeltaR)
    {
      tempParton = partonLHEF;
      partonCounter++;
    }
  }

  for(Candidate *const &parton : partonArray)
  {
    dist = jet->Momentum.DeltaR(parton->Momentum); // take the DR
    isGoodCandidate = true;
    for(Candidate *const &partonLHEF : partonLHEFArray)
    {
      if(parton->Momentum.DeltaR(partonLHEF->Momentum) < 0.01 && parton->PID == partonLHEF->PID && partonLHEF->Charge == parton->Charge)
      {
        isGoodCandidate = false;
        break;
      }
    }

    if(!isGoodCandidate) continue;

    if(parton->D1 != -1 || parton->D2 != -1)
    {
      if((std::abs(parton->PID) < 4 || std::abs(parton->PID) == 21)) continue;
      if(dist < biggerConeSize) contaminations.push_back(parton);
    }
  }

  if(partonCounter != 1)
  {
    jet->FlavorPhys = 0;
  }
  else if(contaminations.empty())
  {
    jet->FlavorPhys = std::abs(tempParton->PID);
  }
  else
  {
    jet->FlavorPhys = std::abs(tempParton->PID);

    for(Candidate *const &parton : contaminations)
    {
      contaminatingFlavor = std::abs(parton->PID);
      motherCounter = 0;
      if(parton->M1 != -1) motherCounter++;
      if(parton->M2 != -1) motherCounter++;

      if(parton->M1 != -1)
      {
        if(const Candidate *mother1 = static_cast<Candidate *>(fParticleInputArray->at(parton->M1));
          mother1 && motherCounter > 0 && mother1->Momentum.DeltaR(tempParton->Momentum) < 0.001) continue;
      }
      if(parton->M2 != -1)
      {
        if(const Candidate *mother2 = static_cast<Candidate *>(fParticleInputArray->at(parton->M2));
          mother2 && motherCounter > 0 && mother2->Momentum.DeltaR(tempParton->Momentum) < 0.001) continue;
      }
      // mother is the initialParton --> OK
      if(std::abs(tempParton->PID) == 4)
      {
        // keep association --> the initialParton is a c --> the contaminated parton is a c
        if(contaminatingFlavor == 4) continue;
        jet->FlavorPhys = 0; // all the other cases reject!
        break;
      }
    }
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("JetFlavorAssociation", JetFlavorAssociation);
