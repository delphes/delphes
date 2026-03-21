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
  EnergySmearing() : fFormula(std::make_unique<DelphesFormula>()) {}

  void Init() override;
  void Process() override;

private:
  const std::unique_ptr<DelphesFormula> fFormula; //!

  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void EnergySmearing::Init()
{
  // read resolution formula
  fFormula->Compile(GetString("ResolutionFormula", "0.0"));

  // import input arrays
  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void EnergySmearing::Process()
{
  fOutputArray->clear();

  Double_t pt, energy, eta, phi, m;

  for(Candidate *const &candidate : *fInputArray)
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;

    pt = candidatePosition.Pt();
    eta = candidatePosition.Eta();
    phi = candidatePosition.Phi();
    energy = candidateMomentum.E();
    m = candidateMomentum.M();

    // apply smearing formula
    energy = gRandom->Gaus(energy, fFormula->Eval(pt, eta, phi, energy));

    if(energy <= 0.0) continue;

    Candidate *new_candidate = static_cast<Candidate *>(candidate->Clone());
    eta = candidateMomentum.Eta();
    phi = candidateMomentum.Phi();
    pt = (energy > m) ? std::sqrt(energy * energy - m * m) / std::cosh(eta) : 0;
    new_candidate->Momentum.SetPtEtaPhiE(pt, eta, phi, energy);
    new_candidate->TrackResolution = fFormula->Eval(pt, eta, phi, energy) / candidateMomentum.E();
    new_candidate->AddCandidate(candidate);

    fOutputArray->emplace_back(new_candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("EnergySmearing", EnergySmearing);
