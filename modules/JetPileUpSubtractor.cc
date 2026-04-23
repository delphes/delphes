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

/** \class JetPileUpSubtractor
 *
 *  Subtract pile-up contribution from jets using the fastjet area method
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>

class JetPileUpSubtractor: public DelphesModule
{
public:
  explicit JetPileUpSubtractor(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fJetPTMin(Steer<double>("JetPTMin", 20.0)) {}

  void Init() override
  {
    fJetInputArray = ImportArray(Steer<std::string>("JetInputArray", "FastJetFinder/jets"));
    fRhoInputArray = ImportArray(Steer<std::string>("RhoInputArray", "Rho/rho"));
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "jets"));
  }
  void Process() override;

private:
  const double fJetPTMin;

  CandidatesCollection fJetInputArray; //!
  CandidatesCollection fRhoInputArray; //!

  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void JetPileUpSubtractor::Process()
{
  fOutputArray->clear();

  // loop over all input candidates
  for(Candidate *const &candidate : *fJetInputArray)
  {
    TLorentzVector momentum = candidate->Momentum;
    const TLorentzVector area = candidate->Area;
    const double eta = momentum.Eta();

    // find rho
    double rho = 0.;
    for(Candidate *const &object : *fRhoInputArray)
    {
      if(eta >= object->Edges[0] && eta < object->Edges[1])
        rho = object->Momentum.Pt();
    }

    // apply pile-up correction
    if(momentum.Pt() <= rho * area.Pt()) continue;

    momentum -= rho * area;

    if(momentum.Pt() <= fJetPTMin) continue;

    Candidate *new_candidate = static_cast<Candidate *>(candidate->Clone());
    new_candidate->Momentum = momentum;

    fOutputArray->emplace_back(new_candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("JetPileUpSubtractor", JetPileUpSubtractor);
