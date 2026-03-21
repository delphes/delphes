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

/** \class MomentumSmearing
 *
 *  Performs transverse momentum resolution smearing.
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

class MomentumSmearing: public DelphesModule
{
public:
  MomentumSmearing() : fFormula(std::make_unique<DelphesFormula>()) {}

  void Init() override;
  void Process() override;

private:
  Double_t LogNormal(Double_t mean, Double_t sigma);

  const std::unique_ptr<DelphesFormula> fFormula; //!

  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!

  Double_t fUseMomentumVector; //!
};

//------------------------------------------------------------------------------

void MomentumSmearing::Init()
{
  // read resolution formula
  fFormula->Compile(GetString("ResolutionFormula", "0.0"));

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));

  // switch to compute momentum smearing based on momentum vector eta, phi
  fUseMomentumVector = GetBool("UseMomentumVector", false);

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void MomentumSmearing::Process()
{
  fOutputArray->clear();

  Double_t pt, eta, phi, e, m, res;

  for(Candidate *const &candidate : *fInputArray)
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    eta = candidatePosition.Eta();
    phi = candidatePosition.Phi();

    if(fUseMomentumVector)
    {
      eta = candidateMomentum.Eta();
      phi = candidateMomentum.Phi();
    }

    pt = candidateMomentum.Pt();
    e = candidateMomentum.E();
    m = candidateMomentum.M();
    res = fFormula->Eval(pt, eta, phi, e, candidate);

    // apply smearing formula
    //pt = gRandom->Gaus(pt, fFormula->Eval(pt, eta, phi, e) * pt);

    res = (res > 1.0) ? 1.0 : res;

    pt = LogNormal(pt, res * pt);

    //if(pt <= 0.0) continue;

    Candidate *new_candidate = static_cast<Candidate *>(candidate->Clone());
    eta = candidateMomentum.Eta();
    phi = candidateMomentum.Phi();
    new_candidate->Momentum.SetPtEtaPhiM(pt, eta, phi, m);
    //new_candidate->TrackResolution = fFormula->Eval(pt, eta, phi, e);
    new_candidate->TrackResolution = res;
    new_candidate->AddCandidate(candidate);

    fOutputArray->emplace_back(new_candidate);
  }
}
//----------------------------------------------------------------

Double_t MomentumSmearing::LogNormal(Double_t mean, Double_t sigma)
{
  if(mean > 0.0)
  {
    const double b = std::sqrt(std::log((1.0 + (sigma * sigma) / (mean * mean)))),
                 a = std::log(mean) - 0.5 * b * b;
    return std::exp(a + b * gRandom->Gaus(0.0, 1.0));
  }
  else
    return 0.0;
}

//------------------------------------------------------------------------------

REGISTER_MODULE("MomentumSmearing", MomentumSmearing);
