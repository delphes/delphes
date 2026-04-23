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
#include <TRandom3.h>

#include "Hector/H_BeamLine.h"
#include "Hector/H_BeamParticle.h"
#include "Hector/H_RecRPObject.h"

using namespace std;

class Hector: public DelphesModule
{
public:
  explicit Hector(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    //
    fDirection(Steer<int>("Direction", 1)),
    fBeamLineLength(Steer<double>("BeamLineLength", 430.0)),
    fDistance(Steer<double>("Distance", 420.0)),
    fOffsetX(Steer<double>("OffsetX", 0.0)),
    fOffsetS(Steer<double>("OffsetS", 120.0)),
    fSigmaE(Steer<double>("SigmaE", 0.0)),
    fSigmaX(Steer<double>("SigmaX", 0.0)),
    fSigmaY(Steer<double>("SigmaY", 0.0)),
    fSigmaT(Steer<double>("SigmaT", 0.0)),
    fEtaMin(Steer<double>("EtaMin", 5.0)),
    fBeamLine(std::make_unique<H_BeamLine>(fDirection, fBeamLineLength + 0.1))
  {
    fBeamLine->fill(Steer<std::string>("BeamLineFile", "cards/LHCB1IR5_5TeV.tfs"), fDirection, moduleParams.Get<std::string>("IPName", "IP5"));
    fBeamLine->offsetElements(fOffsetS, fOffsetX);
    fBeamLine->calcMatrix();
  }

  void Init() override
  {
    fInputArray = ImportArray(Steer<std::string>("InputArray", "ParticlePropagator/stableParticles"));
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "hits"));
  }
  void Process() override;

private:
  const int fDirection;

  const double fBeamLineLength, fDistance;
  const double fOffsetX, fOffsetS;
  const double fSigmaE, fSigmaX, fSigmaY, fSigmaT;
  const double fEtaMin;

  const std::unique_ptr<H_BeamLine> fBeamLine;

  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void Hector::Process()
{
  fOutputArray->clear();

  const double c_light = 2.99792458E8;

  for(Candidate *const &candidate : *fInputArray)
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    const double pz = candidateMomentum.Pz();

    if(std::fabs(candidateMomentum.Eta()) <= fEtaMin || pz * fDirection / std::abs(fDirection) != pz) continue;

    const double x = 1.0E3 * candidatePosition.X(),
                 y = 1.0E3 * candidatePosition.Y(),
                 z = 1.0E-3 * candidatePosition.Z();

    //const double tx = 1.0E6 * std::atan(candidateMomentum.Px() / pz),
    //             ty = 1.0E6 * std::atan(candidateMomentum.Py() / pz);

    const double tx = 0.0, ty = 0.0;

    const double theta = std::hypot(std::atan(candidateMomentum.Px() / pz), std::atan(candidateMomentum.Py() / pz));
    const double distance = (fDistance - 1.0E-3 * candidatePosition.Z()) / std::cos(theta);
    const double time = gRandom->Gaus((distance + 1.0E-3 * candidatePosition.T()) / c_light, fSigmaT);

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

    Candidate *new_candidate = static_cast<Candidate *>(candidate->Clone());
    new_candidate->Position.SetXYZT(particle.getX(), particle.getY(), particle.getS(), time);
    new_candidate->Momentum.SetPxPyPzE(particle.getTX(), particle.getTY(), 0.0, particle.getE());
    new_candidate->AddCandidate(candidate);

    fOutputArray->emplace_back(new_candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("Hector", Hector);
