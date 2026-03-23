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

/** \class PileUpMerger
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

using namespace std;

class PileUpMerger: public DelphesModule
{
public:
  explicit PileUpMerger(const DelphesParameters &moduleParams) :
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
    fReader(std::make_unique<DelphesPileUpReader>(Steer<std::string>("PileUpFile", "MinBias.pileup"))),
    fFunction(std::unique_ptr<DelphesTF2>())
  {
    // read vertex smearing formula
    fFunction->Compile(Steer<std::string>("VertexDistributionFormula", "0.0").data());
    fFunction->SetRange(-fZVertexSpread, -fTVertexSpread, fZVertexSpread, fTVertexSpread);
  }

  void Init() override
  {
    fInputArray = ImportArray(Steer<std::string>("InputArray", "Delphes/stableParticles"));
    fParticleOutputArray = ExportArray(Steer<std::string>("ParticleOutputArray", "stableParticles"));
    fVertexOutputArray = ExportArray(Steer<std::string>("VertexOutputArray", "vertices"));
  }
  void Process() override;

private:
  const int fPileUpDistribution;
  const double fMeanPileUp;

  const double fZVertexSpread;
  const double fTVertexSpread;

  const double fInputBeamSpotX;
  const double fInputBeamSpotY;
  const double fOutputBeamSpotX;
  const double fOutputBeamSpotY;

  const std::unique_ptr<DelphesPileUpReader> fReader; //!
  const std::unique_ptr<DelphesTF2> fFunction; //!

  CandidatesCollection fInputArray; //!

  CandidatesCollection fParticleOutputArray; //!
  CandidatesCollection fVertexOutputArray; //!
};

//------------------------------------------------------------------------------

void PileUpMerger::Process()
{
  fParticleOutputArray->clear();
  fVertexOutputArray->clear();

  TDatabasePDG *pdg = TDatabasePDG::Instance();
  TParticlePDG *pdgParticle;
  int pid, nch, nvtx = -1;
  float x, y, z, t, vx, vy;
  float px, py, pz, e, pt;
  double dz, dphi, dt, sumpt2, dz0, dt0;
  int numberOfEvents, event, numberOfParticles;
  Long64_t allEntries, entry;
  Candidate *candidate, *vertex;
  DelphesFactory *factory;

  const double c_light = 2.99792458E8;

  // --- Deal with primary vertex first  ------

  fFunction->GetRandom2(dz, dt);

  dz0 = -1.0e6;
  dt0 = -1.0e6;

  dt *= c_light * 1.0E3; // necessary in order to make t in mm/c
  dz *= 1.0E3; // necessary in order to make z in mm

  //cout<<dz<<","<<dt<<endl;

  vx = 0.0;
  vy = 0.0;

  numberOfParticles = fInputArray->size();
  nch = 0;
  sumpt2 = 0.0;

  factory = GetFactory();
  vertex = factory->NewCandidate();

  for(Candidate *const &candidate : *fInputArray)
  {
    vx += candidate->Position.X();
    vy += candidate->Position.Y();
    z = candidate->Position.Z();
    t = candidate->Position.T();
    pt = candidate->Momentum.Pt();

    // take postion and time from first stable particle
    if(dz0 < -999999.0)
      dz0 = z;
    if(dt0 < -999999.0)
      dt0 = t;

    // cancel any possible offset in position and time the input file
    candidate->Position.SetZ(z - dz0 + dz);
    candidate->Position.SetT(t - dt0 + dt);

    candidate->IsPU = 0;

    fParticleOutputArray->emplace_back(candidate);

    if(std::fabs(candidate->Charge) > 1.0E-9)
    {
      nch++;
      sumpt2 += pt * pt;
      vertex->AddCandidate(candidate);
    }
  }

  if(numberOfParticles > 0)
  {
    vx /= sumpt2;
    vy /= sumpt2;
  }

  nvtx++;
  vertex->Position.SetXYZT(vx, vy, dz, dt);
  vertex->ClusterIndex = nvtx;
  vertex->ClusterNDF = nch;
  vertex->SumPT2 = sumpt2;
  vertex->GenSumPT2 = sumpt2;
  fVertexOutputArray->emplace_back(vertex);

  // --- Then with pile-up vertices  ------

  switch(fPileUpDistribution)
  {
  case 0:
    numberOfEvents = gRandom->Poisson(fMeanPileUp);
    break;
  case 1:
    numberOfEvents = gRandom->Integer(2 * fMeanPileUp + 1);
    break;
  case 2:
    numberOfEvents = fMeanPileUp;
    break;
  default:
    numberOfEvents = gRandom->Poisson(fMeanPileUp);
    break;
  }

  allEntries = fReader->GetEntries();

  for(event = 0; event < numberOfEvents; ++event)
  {
    do
    {
      entry = std::ceil(gRandom->Rndm() * allEntries);
    } while(entry >= allEntries);

    fReader->ReadEntry(entry);

    // --- Pile-up vertex smearing

    fFunction->GetRandom2(dz, dt);

    dt *= c_light * 1.0E3; // necessary in order to make t in mm/c
    dz *= 1.0E3; // necessary in order to make z in mm

    dphi = gRandom->Uniform(-M_PI, M_PI);

    vx = 0.0;
    vy = 0.0;

    numberOfParticles = 0;
    sumpt2 = 0.0;

    //factory = GetFactory();
    vertex = factory->NewCandidate();

    while(fReader->ReadParticle(pid, x, y, z, t, px, py, pz, e))
    {
      candidate = factory->NewCandidate();

      candidate->PID = pid;

      candidate->Status = 1;

      pdgParticle = pdg->GetParticle(pid);
      candidate->Charge = pdgParticle ? int(pdgParticle->Charge() / 3.0) : -999;
      candidate->Mass = pdgParticle ? pdgParticle->Mass() : -999.9;

      candidate->IsPU = 1;

      candidate->Momentum.SetPxPyPzE(px, py, pz, e);
      candidate->Momentum.RotateZ(dphi);
      pt = candidate->Momentum.Pt();

      x -= fInputBeamSpotX;
      y -= fInputBeamSpotY;
      candidate->Position.SetXYZT(x, y, z + dz, t + dt);
      candidate->Position.RotateZ(dphi);
      candidate->Position += TLorentzVector(fOutputBeamSpotX, fOutputBeamSpotY, 0.0, 0.0);

      vx += candidate->Position.X();
      vy += candidate->Position.Y();

      ++numberOfParticles;
      if(std::fabs(candidate->Charge) > 1.0E-9)
      {
        nch++;
        sumpt2 += pt * pt;
        vertex->AddCandidate(candidate);
      }

      fParticleOutputArray->emplace_back(candidate);
    }

    if(numberOfParticles > 0)
    {
      vx /= numberOfParticles;
      vy /= numberOfParticles;
    }

    nvtx++;

    vertex->Position.SetXYZT(vx, vy, dz, dt);

    vertex->ClusterIndex = nvtx;
    vertex->ClusterNDF = nch;
    vertex->SumPT2 = sumpt2;
    vertex->GenSumPT2 = sumpt2;

    vertex->IsPU = 1;

    fVertexOutputArray->emplace_back(vertex);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("PileUpMerger", PileUpMerger);
