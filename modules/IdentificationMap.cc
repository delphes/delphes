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

/** \class IdentificationMap
 *
 *  Converts particles with some PDG code into another particle,
 *  according to parametrized probability.
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>
#include <TRandom3.h>

using namespace std;

class IdentificationMap: public DelphesModule
{
public:
  explicit IdentificationMap(const DelphesParameters &moduleParams) : DelphesModule(moduleParams)
  {
    // read efficiency formulas
    for(const std::pair<int, std::pair<int, std::string> > efficiencyFormula :
      Steer<std::unordered_map<int, std::pair<int, std::string> > >("EfficiencyFormula"))
    {
      std::unique_ptr<DelphesFormula> formula = std::make_unique<DelphesFormula>();
      formula->Compile(efficiencyFormula.second.second);
      const int pdg = efficiencyFormula.first;
      const int pdg2 = efficiencyFormula.second.first;
      fEfficiencyMap.insert(std::make_pair(pdg, make_pair(pdg2, std::move(formula))));
    }
    // set default efficiency formula
    if(fEfficiencyMap.count(0) == 0)
    {
      std::unique_ptr<DelphesFormula> formula = std::make_unique<DelphesFormula>();
      formula->Compile("1.0");
      fEfficiencyMap.insert(std::make_pair(0, std::make_pair(0, std::move(formula))));
    }
  }

  void Init() override
  {
    fInputArray = ImportArray(Steer<std::string>("InputArray", "ParticlePropagator/stableParticles"));
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "stableParticles"));
  }
  void Process() override;

private:
  typedef std::multimap<Int_t, std::pair<Int_t, std::unique_ptr<DelphesFormula> > > TMisIDMap; //!

  TMisIDMap fEfficiencyMap; //!

  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void IdentificationMap::Process()
{
  fOutputArray->clear();

  for(Candidate *const &candidate : *fInputArray)
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    const double eta = candidatePosition.Eta();
    const double phi = candidatePosition.Phi();
    const double pt = candidateMomentum.Pt();
    const double e = candidateMomentum.E();

    const int pdgCodeIn = candidate->PID;
    const int charge = candidate->Charge;

    // first check that PID of this particle is specified in the map
    // otherwise, look for PID = 0

    //TMisIDMap::iterator itEfficiencyMap = fEfficiencyMap.find(pdgCodeIn);

    std::pair<TMisIDMap::iterator, TMisIDMap::iterator> range = fEfficiencyMap.equal_range(pdgCodeIn);
    if(range.first == range.second) range = fEfficiencyMap.equal_range(-pdgCodeIn);
    if(range.first == range.second) range = fEfficiencyMap.equal_range(0);

    const double r = gRandom->Uniform();
    double total = 0.;

    // loop over sub-map for this PID
    for(TMisIDMap::iterator it = range.first; it != range.second; ++it)
    {
      std::unique_ptr<DelphesFormula> &formula = (it->second).second;
      const int pdgCodeOut = (it->second).first;

      const double p = formula->Eval(pt, eta, phi, e);

      if(total <= r && r < total + p)
      {
        // change PID of particle
        Candidate *new_candidate = static_cast<Candidate *>(candidate->Clone());
        if(pdgCodeOut != 0) new_candidate->PID = charge * pdgCodeOut;
        fOutputArray->emplace_back(new_candidate);
        break;
      }

      total += p;
    }
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("IdentificationMap", IdentificationMap);
