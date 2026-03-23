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

/** \class PileUpMergerPythia8
 *
 *  Merges particles from pile-up sample into event
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesModule.h"
#include "classes/DelphesPileUpReader.h"
#include "classes/DelphesTF2.h"

#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TRandom3.h>

#include <Pythia.h>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

class PileUpMergerPythia8: public DelphesModule
{
public:
  explicit PileUpMergerPythia8(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fPileUpDistribution(Steer<int>("PileUpDistribution", 0)),
    fMeanPileUp(Steer<double>("MeanPileUp", 10)),
    //
    fZVertexSpread(Steer<double>("ZVertexSpread", 0.15)),
    fTVertexSpread(Steer<double>("TVertexSpread", 1.5E-09)),
    //
    fInputBeamSpotX(Steer<double>("InputBeamSpotX", 0.0)),
    fInputBeamSpotY(Steer<double>("InputBeamSpotY", 0.0)),
    fOutputBeamSpotX(Steer<double>("OutputBeamSpotX", 0.0)),
    fOutputBeamSpotY(Steer<double>("OutputBeamSpotY", 0.0)),
    //
    fPTMin(Steer<double>("PTMin", 0.0)),
    //
    fConfigFile(Steer<std::string>("ConfigFile", "MinBias.cmnd")),
    //
    fFunction(std::make_unique<DelphesTF2>())
  {
    fFunction->Compile(Steer<std::string>("VertexDistributionFormula", "0.0").data());
    fFunction->SetRange(-fZVertexSpread, -fTVertexSpread, fZVertexSpread, fTVertexSpread);
  }

  void Init() override
  {
    fPythia = std::make_unique<Pythia8::Pythia>();
    fPythia->readFile(fConfigFile);

    fInputArray = ImportArray(Steer<std::string>("InputArray", "Delphes/stableParticles"));
    fParticleOutputArray = ExportArray(Steer<std::string>("ParticleOutputArray", "stableParticles"));
    fVertexOutputArray = ExportArray(Steer<std::string>("VertexOutputArray", "vertices"));
  }
  void Process() override;
  void Finish() override
  {
    fPythia.reset();
  }

private:
  const int fPileUpDistribution;
  const double fMeanPileUp;

  const double fZVertexSpread;
  const double fTVertexSpread;

  const double fInputBeamSpotX;
  const double fInputBeamSpotY;
  const double fOutputBeamSpotX;
  const double fOutputBeamSpotY;

  const double fPTMin;

  const std::string fConfigFile;

  const std::unique_ptr<DelphesTF2> fFunction; //!

  std::unique_ptr<Pythia8::Pythia> fPythia; //!

  CandidatesCollection fInputArray; //!

  CandidatesCollection fParticleOutputArray; //!
  CandidatesCollection fVertexOutputArray; //!
};

//------------------------------------------------------------------------------

void PileUpMergerPythia8::Process()
{
  fParticleOutputArray->clear();
  fVertexOutputArray->clear();

  TDatabasePDG *pdg = TDatabasePDG::Instance();
  DelphesFactory *factory = GetFactory();

  const double c_light = 2.99792458E8;

  // --- Deal with primary vertex first  ------

  double dz, dt;
  fFunction->GetRandom2(dz, dt);

  dt *= c_light * 1.0E3; // necessary in order to make t in mm/c
  dz *= 1.0E3; // necessary in order to make z in mm
  float vx = 0.f, vy = 0.f;
  const size_t numberOfParticles = fInputArray->size();
  for(Candidate *const &candidate : *fInputArray)
  {
    vx += candidate->Position.X();
    vy += candidate->Position.Y();
    const float z = candidate->Position.Z(), t = candidate->Position.T();
    candidate->Position.SetZ(z + dz);
    candidate->Position.SetT(t + dt);
    fParticleOutputArray->emplace_back(candidate);
  }

  if(numberOfParticles > 0)
  {
    vx /= numberOfParticles;
    vy /= numberOfParticles;
  }

  Candidate *vertex = factory->NewCandidate();
  vertex->Position.SetXYZT(vx, vy, dz, dt);
  fVertexOutputArray->emplace_back(vertex);

  // --- Then with pile-up vertices  ------

  int numberOfEvents = 0;
  switch(fPileUpDistribution)
  {
  case 0:
    numberOfEvents = gRandom->Poisson(fMeanPileUp);
    break;
  case 1:
    numberOfEvents = gRandom->Integer(2 * fMeanPileUp + 1);
    break;
  default:
    numberOfEvents = gRandom->Poisson(fMeanPileUp);
    break;
  }

  for(int event = 0; event < numberOfEvents; ++event)
  {
    while(!fPythia->next());

    // --- Pile-up vertex smearing

    fFunction->GetRandom2(dz, dt);

    dt *= c_light * 1.0E3; // necessary in order to make t in mm/c
    dz *= 1.0E3; // necessary in order to make z in mm

    const double dphi = gRandom->Uniform(-M_PI, M_PI);

    float vx = 0.f, vy = 0.f;
    const size_t numberOfParticles = fPythia->event.size();
    for(size_t i = 1; i < numberOfParticles; ++i)
    {
      Pythia8::Particle &particle = fPythia->event[i];

      const int status = particle.statusHepMC();

      if(status != 1 || !particle.isVisible() || particle.pT() <= fPTMin) continue;

      const int pid = particle.id();
      const float px = particle.px(), py = particle.py(), pz = particle.pz(), e = particle.e();
      float x = particle.xProd(), y = particle.yProd(), z = particle.zProd(), t = particle.tProd();

      Candidate *candidate = factory->NewCandidate();

      candidate->PID = pid;

      candidate->Status = 1;

      TParticlePDG *pdgParticle = pdg->GetParticle(pid);
      candidate->Charge = pdgParticle ? int(pdgParticle->Charge() / 3.0) : -999;
      candidate->Mass = pdgParticle ? pdgParticle->Mass() : -999.9;

      candidate->IsPU = 1;

      candidate->Momentum.SetPxPyPzE(px, py, pz, e);
      candidate->Momentum.RotateZ(dphi);

      x -= fInputBeamSpotX;
      y -= fInputBeamSpotY;
      candidate->Position.SetXYZT(x, y, z + dz, t + dt);
      candidate->Position.RotateZ(dphi);
      candidate->Position += TLorentzVector(fOutputBeamSpotX, fOutputBeamSpotY, 0.0, 0.0);

      vx += candidate->Position.X();
      vy += candidate->Position.Y();

      fParticleOutputArray->emplace_back(candidate);
    }

    if(numberOfParticles > 0)
    {
      vx /= numberOfParticles;
      vy /= numberOfParticles;
    }

    Candidate *vertex = factory->NewCandidate();
    vertex->Position.SetXYZT(vx, vy, dz, dt);
    vertex->IsPU = 1;

    fVertexOutputArray->emplace_back(vertex);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("PileUpMergerPythia8", PileUpMergerPythia8);
