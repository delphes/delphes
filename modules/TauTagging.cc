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

#include "modules/TauTaggingPartonClassifier.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFilter.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>
#include <TRandom3.h>

#include <map>

using namespace std;

class TauTagging: public DelphesModule
{
public:
  explicit TauTagging(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fBitNumber(Steer<int>("BitNumber", 0)),
    fDeltaR(Steer<double>("DeltaR", 0.5))
  {
    for(const std::pair<int, std::string> &efficiencyFormula :
      Steer<std::vector<std::pair<int, std::string> > >("EfficiencyFormula"))
    {
      std::unique_ptr<DelphesFormula> formula = std::make_unique<DelphesFormula>();
      formula->Compile(efficiencyFormula.second);
      fEfficiencyMap[efficiencyFormula.first] = std::move(formula);
    }
    // set default efficiency formula
    if(fEfficiencyMap.count(0) == 0)
    {
      std::unique_ptr<DelphesFormula> formula = std::make_unique<DelphesFormula>();
      formula->Compile("0.0");
      fEfficiencyMap[0] = std::move(formula);
    }
  }

  void Init() override
  {
    fParticleInputArray = ImportArray(Steer<std::string>("ParticleInputArray", "Delphes/allParticles"));
    fPartonInputArray = ImportArray(Steer<std::string>("PartonInputArray", "Delphes/partons"));
    fJetInputArray = ImportArray(Steer<std::string>("JetInputArray", "FastJetFinder/jets"));

    fClassifier = std::make_unique<TauTaggingPartonClassifier>(fParticleInputArray);
    fClassifier->fPTMin = Steer<double>("TauPTMin", 1.0);
    fClassifier->fEtaMax = Steer<double>("TauEtaMax", 2.5);

    fFilter = std::make_unique<DelphesFilter>(fPartonInputArray);
  }
  void Process() override;

private:
  const int fBitNumber;
  const double fDeltaR;

  CandidatesCollection fParticleInputArray; //!
  CandidatesCollection fPartonInputArray; //!
  CandidatesCollection fJetInputArray; //!

  std::unique_ptr<TauTaggingPartonClassifier> fClassifier; //!
  std::unique_ptr<DelphesFilter> fFilter;

#if !defined(__CINT__) && !defined(__CLING__)
  std::map<int, std::unique_ptr<DelphesFormula> > fEfficiencyMap; //!
#endif
};

//------------------------------------------------------------------------------

void TauTagging::Process()
{
  // select taus
  fFilter->Reset();
  const std::vector<Candidate *> tauArray = fFilter->GetSubArray(fClassifier.get(), 0);

  // loop over all input jets
  for(Candidate *const &jet : *fJetInputArray)
  {
    const TLorentzVector &jetMomentum = jet->Momentum;
    int pdgCode = 0;
    int charge = gRandom->Uniform() > 0.5 ? 1 : -1;
    const double eta = jetMomentum.Eta();
    const double phi = jetMomentum.Phi();
    const double pt = jetMomentum.Pt();
    const double e = jetMomentum.E();

    // loop over all input taus
    for(Candidate *const &tau : tauArray)
    {
      if(tau->D1 < 0) continue;

      if(tau->D1 >= static_cast<int>(fParticleInputArray->size()) || tau->D2 >= static_cast<int>(fParticleInputArray->size()))
      {
        throw runtime_error("tau's daughter index is greater than the ParticleInputArray size");
      }

      TLorentzVector tauMomentum;
      for(int i = tau->D1; i <= tau->D2; ++i)
      {
        Candidate *const &daughter = static_cast<Candidate *>(fParticleInputArray->at(i));
        if(std::abs(daughter->PID) == 16) continue;
        tauMomentum += daughter->Momentum;
      }

      if(jetMomentum.DeltaR(tauMomentum) <= fDeltaR)
      {
        pdgCode = 15;
        charge = tau->Charge;
      }
    }

    // fake electrons and muons

    if(pdgCode == 0)
    {

      double drMin = fDeltaR;
      for(Candidate *const &part : *fPartonInputArray)
      {
        if(std::abs(part->PID) == 11 || std::abs(part->PID) == 13)
        {
          TLorentzVector &tauMomentum = part->Momentum;
          if(tauMomentum.Pt() < fClassifier->fPTMin) continue;
          if(std::fabs(tauMomentum.Eta()) > fClassifier->fEtaMax) continue;

          double dr = jetMomentum.DeltaR(tauMomentum);
          if(dr < drMin)
          {
            drMin = dr;
            pdgCode = std::abs(part->PID);
            charge = part->Charge;
          }
        }
      }
    }

    // find an efficency formula
    std::unique_ptr<DelphesFormula> &formula =
      fEfficiencyMap.count(pdgCode) == 0 ? fEfficiencyMap.at(0) : fEfficiencyMap.at(pdgCode);

    // apply an efficency formula
    const double eff = formula->Eval(pt, eta, phi, e);
    jet->TauFlavor = pdgCode;
    jet->TauTag |= (gRandom->Uniform() <= eff) << fBitNumber;
    jet->TauWeight = eff;

    // set tau charge
    jet->Charge = charge;
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("TauTagging", TauTagging);
