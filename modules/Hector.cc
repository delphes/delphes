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

/** \class Hector
 *
 *  Propagates candidates using Hector library.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>

#include "Hector/H_BeamLine.h"
#include "Hector/H_BeamParticle.h"
#include "Hector/H_RecRPObject.h"

using namespace std;

class Hector: public DelphesModule
{
public:
  Hector() = default;

  void Init() override;
  void Process() override;

private:
  Int_t fDirection;

  Double_t fBeamLineLength, fDistance;
  Double_t fOffsetX, fOffsetS;
  Double_t fSigmaE, fSigmaX, fSigmaY, fSigmaT;
  Double_t fEtaMin;

  std::unique_ptr<H_BeamLine> fBeamLine;

  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void Hector::Init()
{
  // read Hector parameters

  fDirection = GetInt("Direction", 1);
  fBeamLineLength = GetDouble("BeamLineLength", 430.0);
  fDistance = GetDouble("Distance", 420.0);
  fOffsetX = GetDouble("OffsetX", 0.0);
  fOffsetS = GetDouble("OffsetS", 120.0);
  fSigmaE = GetDouble("SigmaE", 0.0);
  fSigmaX = GetDouble("SigmaX", 0.0);
  fSigmaY = GetDouble("SigmaY", 0.0);
  fSigmaT = GetDouble("SigmaT", 0.0);
  fEtaMin = GetDouble("EtaMin", 5.0);

  fBeamLine = std::make_unique<H_BeamLine>(fDirection, fBeamLineLength + 0.1);
  fBeamLine->fill(GetString("BeamLineFile", "cards/LHCB1IR5_5TeV.tfs"), fDirection, GetString("IPName", "IP5"));
  fBeamLine->offsetElements(fOffsetS, fOffsetX);
  fBeamLine->calcMatrix();

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "hits"));
}

//------------------------------------------------------------------------------

void Hector::Process()
{
  fOutputArray->clear();

  Double_t pz;
  Double_t x, y, z, tx, ty, theta;
  Double_t distance, time;

  const Double_t c_light = 2.99792458E8;

  for(const auto &candidate : *fInputArray)
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    pz = candidateMomentum.Pz();

    if(TMath::Abs(candidateMomentum.Eta()) <= fEtaMin || TMath::Sign(pz, Double_t(fDirection)) != pz) continue;

    x = 1.0E3 * candidatePosition.X();
    y = 1.0E3 * candidatePosition.Y();
    z = 1.0E-3 * candidatePosition.Z();

    //    tx = 1.0E6 * TMath::ATan(candidateMomentum.Px()/pz);
    //    ty = 1.0E6 * TMath::ATan(candidateMomentum.Py()/pz);

    tx = 0.0;
    ty = 0.0;

    theta = TMath::Hypot(TMath::ATan(candidateMomentum.Px() / pz), TMath::ATan(candidateMomentum.Py() / pz));
    distance = (fDistance - 1.0E-3 * candidatePosition.Z()) / TMath::Cos(theta);
    time = gRandom->Gaus((distance + 1.0E-3 * candidatePosition.T()) / c_light, fSigmaT);

    H_BeamParticle particle(candidate->Mass, candidate->Charge);
    //    particle.set4Momentum(candidateMomentum);
    particle.set4Momentum(candidateMomentum.Px(), candidateMomentum.Py(),
      candidateMomentum.Pz(), candidateMomentum.E());
    particle.setPosition(x, y, tx, ty, z);

    particle.smearAng(fSigmaX, fSigmaY, gRandom);
    particle.smearE(fSigmaE, gRandom);

    particle.computePath(fBeamLine.get());

    if(particle.stopped(fBeamLine.get())) continue;

    particle.propagate(fDistance);

    auto *new_candidate = static_cast<Candidate *>(candidate->Clone());
    new_candidate->Position.SetXYZT(particle.getX(), particle.getY(), particle.getS(), time);
    new_candidate->Momentum.SetPxPyPzE(particle.getTX(), particle.getTY(), 0.0, particle.getE());
    new_candidate->AddCandidate(candidate);

    fOutputArray->emplace_back(new_candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("Hector", Hector);
