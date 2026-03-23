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

/** \class EnergySmearing
 *
 *  Performs energy resolution smearing.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>
#include <TRandom3.h>

using namespace std;

class EnergySmearing: public DelphesModule
{
public:
  explicit EnergySmearing(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fFormula(std::make_unique<DelphesFormula>())
  {
    // read resolution formula
    fFormula->Compile(Steer<std::string>("ResolutionFormula", "0.0"));
  }

  void Init() override
  {
    fInputArray = ImportArray(Steer<std::string>("InputArray", "ParticlePropagator/stableParticles")); // import input arrays
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "stableParticles")); // create output array
  }
  void Process() override;

private:
  const std::unique_ptr<DelphesFormula> fFormula; //!

  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void EnergySmearing::Process()
{
  fOutputArray->clear();

  for(Candidate *const &candidate : *fInputArray)
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;

    double energy = candidateMomentum.E();

    // apply smearing formula
    energy = gRandom->Gaus(energy,
      fFormula->Eval(candidatePosition.Pt(), candidatePosition.Eta(), candidatePosition.Phi(), energy));

    if(energy <= 0.0) continue;

    Candidate *new_candidate = static_cast<Candidate *>(candidate->Clone());
    const double eta = candidateMomentum.Eta();
    const double phi = candidateMomentum.Phi();
    const double m = candidateMomentum.M();
    const double pt = (energy > m) ? std::sqrt(energy * energy - m * m) / std::cosh(eta) : 0;
    new_candidate->Momentum.SetPtEtaPhiE(pt, eta, phi, energy);
    new_candidate->TrackResolution = fFormula->Eval(pt, eta, phi, energy) / candidateMomentum.E();
    new_candidate->AddCandidate(candidate);

    fOutputArray->emplace_back(new_candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("EnergySmearing", EnergySmearing);
