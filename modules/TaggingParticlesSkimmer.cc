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

/** \class TaggingParticlesSkimmer
 *
 *  Filters particle collection by only keeping gen particles useful for b and tau tagging.
    tau particles are replaced by their "visible decay".
 *
 *  \author M. Selvaggi
 *
 */

#include "modules/TauTaggingPartonClassifier.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFilter.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>

using namespace std;

class TaggingParticlesSkimmer: public DelphesModule
{
public:
  explicit TaggingParticlesSkimmer(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fPTMin(Steer<double>("PTMin", 15.0)),
    fEtaMax(Steer<double>("EtaMax", 2.5)) {}

  void Init() override
  {
    fPartonInputArray = ImportArray(Steer<std::string>("PartonInputArray", "Delphes/partons"));
    fParticleInputArray = ImportArray(Steer<std::string>("ParticleInputArray", "Delphes/allParticles"));
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "taggingParticles"));

    fClassifier = std::make_unique<TauTaggingPartonClassifier>(fParticleInputArray);
    fClassifier->fPTMin = Steer<double>("PTMin", 15.0);
    fClassifier->fEtaMax = Steer<double>("EtaMax", 2.5);

    fFilter = std::make_unique<DelphesFilter>(fPartonInputArray);
  }
  void Process() override;

private:
  const double fPTMin; //!
  const double fEtaMax; //!

  CandidatesCollection fPartonInputArray; //!
  CandidatesCollection fParticleInputArray; //!

  CandidatesCollection fOutputArray; //!

  std::unique_ptr<TauTaggingPartonClassifier> fClassifier; //!
  std::unique_ptr<DelphesFilter> fFilter;
};

//------------------------------------------------------------------------------

void TaggingParticlesSkimmer::Process()
{
  fOutputArray->clear();

  // first select hadronic taus and replace them by visible part
  fFilter->Reset();
  const CandidatesCollection tauArray = fFilter->GetSubArray(fClassifier.get(), 0);

  if(tauArray->empty()) return;

  // loop over all input taus
  for(Candidate *const &tau : *tauArray)
  {
    if(tau->D1 < 0) continue;

    if(tau->D1 >= static_cast<int>(fParticleInputArray->size()) || tau->D2 >= static_cast<int>(fParticleInputArray->size()))
    {
      throw runtime_error("tau's daughter index is greater than the ParticleInputArray size");
    }

    TLorentzVector tauMomentum;
    for(int i = tau->D1; i <= tau->D2; ++i)
    {
      Candidate *daughter = static_cast<Candidate *>(fParticleInputArray->at(i));
      if(std::abs(daughter->PID) == 16) continue;
      tauMomentum += daughter->Momentum;
    }

    Candidate *candidate = static_cast<Candidate *>(tau->Clone());
    candidate->Momentum = tauMomentum;

    fOutputArray->emplace_back(candidate);
  }

  // then add all other partons (except tau's to avoid double counting)

  for(Candidate *const &candidate : *fPartonInputArray)
  {
    const int pdgCode = std::abs(candidate->PID);
    if(pdgCode == 15) continue;

    const double pt = candidate->Momentum.Pt();
    if(pt < fPTMin) continue;

    const double eta = std::fabs(candidate->Momentum.Eta());
    if(eta > fEtaMax) continue;

    fOutputArray->emplace_back(candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("TaggingParticlesSkimmer", TaggingParticlesSkimmer);
