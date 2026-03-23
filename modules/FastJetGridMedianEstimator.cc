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

/** \class FastJetGridMedianEstimator
 *
 *  Computes median energy density per event using a fixed grid.
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>

#include <fastjet/tools/GridMedianBackgroundEstimator.hh>

class FastJetGridMedianEstimator: public DelphesModule
{
public:
  explicit FastJetGridMedianEstimator(const DelphesParameters &moduleParams) : DelphesModule(moduleParams)
  {
    for(const std::array<double, 4> &gridRange : Steer<std::vector<std::array<double, 4> > >("GridRange"))
    {
      const double rapMin = gridRange.at(0), rapMax = gridRange.at(1), drap = gridRange.at(2), dphi = gridRange.at(3);
      fEstimators.push_back(std::make_unique<fastjet::GridMedianBackgroundEstimator>(rapMin, rapMax, drap, dphi));
    }
  }

  void Init() override
  {
    fInputArray = ImportArray(Steer<std::string>("InputArray", "Calorimeter/towers"));
    fRhoOutputArray = ExportArray(Steer<std::string>("RhoOutputArray", "rho"));
  }
  void Process() override;

private:
  std::vector<std::unique_ptr<fastjet::GridMedianBackgroundEstimator> > fEstimators; //!

  CandidatesCollection fInputArray;
  CandidatesCollection fRhoOutputArray; //!
};

//------------------------------------------------------------------------------

void FastJetGridMedianEstimator::Process()
{
  fRhoOutputArray->clear();

  DelphesFactory *factory = GetFactory();

  std::vector<fastjet::PseudoJet> inputList;
  size_t number = 0;
  for(Candidate *const &candidate : *fInputArray) // loop over input objects
  {
    TLorentzVector &momentum = candidate->Momentum;
    fastjet::PseudoJet &jet = inputList.emplace_back(momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E());
    jet.set_user_index(number++);
  }

  // compute rho and store it
  for(std::vector<std::unique_ptr<fastjet::GridMedianBackgroundEstimator> >::iterator itEstimators = fEstimators.begin();
    itEstimators != fEstimators.end(); ++itEstimators)
  {
    (*itEstimators)->set_particles(inputList);

    double rho = (*itEstimators)->rho();

    Candidate *candidate = factory->NewCandidate();
    candidate->Momentum.SetPtEtaPhiE(rho, 0.0, 0.0, rho);
    candidate->Edges[0] = (*itEstimators)->rapmin();
    candidate->Edges[1] = (*itEstimators)->rapmax();
    fRhoOutputArray->emplace_back(candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("FastJetGridMedianEstimator", FastJetGridMedianEstimator);
