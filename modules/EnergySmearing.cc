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
#include "classes/DelphesModuleFactory.h"

#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
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

  const TObjArray *fInputArray{nullptr}; //!
  std::unique_ptr<TIterator> fItInputArray; //!

  TObjArray *fOutputArray{nullptr}; //!
};

//------------------------------------------------------------------------------

void EnergySmearing::Init()
{
  // read resolution formula
  fFormula->Compile(GetString("ResolutionFormula", "0.0"));

  // import input arrays
  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));
  fItInputArray.reset(fInputArray->MakeIterator());

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void EnergySmearing::Process()
{
  Candidate *candidate = nullptr, *mother = nullptr;
  Double_t pt, energy, eta, phi, m;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
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

    mother = candidate;
    candidate = static_cast<Candidate *>(candidate->Clone());
    eta = candidateMomentum.Eta();
    phi = candidateMomentum.Phi();
    pt = (energy > m) ? TMath::Sqrt(energy * energy - m * m) / TMath::CosH(eta) : 0;
    candidate->Momentum.SetPtEtaPhiE(pt, eta, phi, energy);
    candidate->TrackResolution = fFormula->Eval(pt, eta, phi, energy) / candidateMomentum.E();
    candidate->AddCandidate(mother);

    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("EnergySmearing", EnergySmearing);
