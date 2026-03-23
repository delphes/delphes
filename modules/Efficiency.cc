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

/** \class Efficiency
 *
 *  Selects candidates from the InputArray according to the efficiency formula.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>
#include <TRandom3.h>

using namespace std;

class Efficiency: public DelphesModule
{
public:
  explicit Efficiency(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fUseMomentumVector(Steer<bool>("UseMomentumVector", false)), // switch to compute efficiency based on momentum vector eta, phi
    fFormula(std::make_unique<DelphesFormula>())
  {
    // read efficiency formula
    fFormula->Compile(Steer<std::string>("EfficiencyFormula", "1.0"));
  }

  void Init() override
  {
    fInputArray = ImportArray(Steer<std::string>("InputArray", "ParticlePropagator/stableParticles")); // import input array
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "stableParticles")); // create output array
  }
  void Process() override;

private:
  const Double_t fUseMomentumVector; //!

  const std::unique_ptr<DelphesFormula> fFormula; //!

  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void Efficiency::Process()
{
  fOutputArray->clear();
  for(Candidate *const &candidate : *fInputArray)
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    double eta = candidatePosition.Eta();
    double phi = candidatePosition.Phi();

    if(fUseMomentumVector)
    {
      eta = candidateMomentum.Eta();
      phi = candidateMomentum.Phi();
    }

    const double pt = candidateMomentum.Pt();
    const double e = candidateMomentum.E();

    // apply an efficency formula
    if(gRandom->Uniform() > fFormula->Eval(pt, eta, phi, e, candidate)) continue;

    fOutputArray->emplace_back(candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("Efficiency", Efficiency);
