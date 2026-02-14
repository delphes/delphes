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

#include "modules/PileUpMergerPythia8.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesPileUpReader.h"
#include "classes/DelphesTF2.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "Pythia.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TString.h"
#include <Math/RotationZ.h>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

PileUpMergerPythia8::PileUpMergerPythia8() :
  fFunction(0), fPythia(0)
{
  fFunction = new DelphesTF2;
}

//------------------------------------------------------------------------------

PileUpMergerPythia8::~PileUpMergerPythia8()
{
  delete fFunction;
}

//------------------------------------------------------------------------------

void PileUpMergerPythia8::Init()
{
  const char *fileName;

  fPileUpDistribution = GetInt("PileUpDistribution", 0);

  fMeanPileUp = GetDouble("MeanPileUp", 10);

  fZVertexSpread = GetDouble("ZVertexSpread", 0.15);
  fTVertexSpread = GetDouble("TVertexSpread", 1.5E-09);

  fInputBeamSpotX = GetDouble("InputBeamSpotX", 0.0);
  fInputBeamSpotY = GetDouble("InputBeamSpotY", 0.0);
  fOutputBeamSpotX = GetDouble("OutputBeamSpotX", 0.0);
  fOutputBeamSpotY = GetDouble("OutputBeamSpotY", 0.0);

  fPTMin = GetDouble("PTMin", 0.0);

  fFunction->Compile(GetString("VertexDistributionFormula", "0.0"));
  fFunction->SetRange(-fZVertexSpread, -fTVertexSpread, fZVertexSpread, fTVertexSpread);

  fileName = GetString("ConfigFile", "MinBias.cmnd");
  fPythia = new Pythia8::Pythia();
  fPythia->readFile(fileName);

  // import input array
  GetFactory()->EventModel()->Attach(GetString("InputArray", "Delphes/stableParticles"), fInputArray); // I/O

  // create output arrays
  GetFactory()->EventModel()->Book(fParticleOutputArray, GetString("ParticleOutputArray", "stableParticles"));
  GetFactory()->EventModel()->Book(fVertexOutputArray, GetString("VertexOutputArray", "vertices"));
}

//------------------------------------------------------------------------------

void PileUpMergerPythia8::Finish()
{
  if(fPythia) delete fPythia;
}

//------------------------------------------------------------------------------

void PileUpMergerPythia8::Process()
{
  TDatabasePDG *pdg = TDatabasePDG::Instance();
  TParticlePDG *pdgParticle;
  Int_t pid, status;
  Float_t x, y, z, t, vx, vy;
  Float_t px, py, pz, e;
  Double_t dz, dphi, dt;
  Int_t numberOfEvents, event, numberOfParticles, i;
  DelphesFactory *factory;

  const Double_t c_light = 2.99792458E8;

  // --- Deal with primary vertex first  ------

  fFunction->GetRandom2(dz, dt);

  dt *= c_light * 1.0E3; // necessary in order to make t in mm/c
  dz *= 1.0E3; // necessary in order to make z in mm
  vx = 0.0;
  vy = 0.0;
  numberOfParticles = fInputArray->size();
  for(auto &candidate : *fInputArray)
  {
    vx += candidate.Position.X();
    vy += candidate.Position.Y();
    z = candidate.Position.Z();
    t = candidate.Position.T();
    candidate.Position.SetPz(z + dz);
    candidate.Position.SetE(t + dt);
    fParticleOutputArray->emplace_back(candidate);
  }

  if(numberOfParticles > 0)
  {
    vx /= numberOfParticles;
    vy /= numberOfParticles;
  }

  factory = GetFactory();

  auto *vertex = factory->NewCandidate();
  vertex->Position.SetXYZT(vx, vy, dz, dt);
  fVertexOutputArray->emplace_back(*vertex);

  // --- Then with pile-up vertices  ------

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

  for(event = 0; event < numberOfEvents; ++event)
  {
    while(!fPythia->next());

    // --- Pile-up vertex smearing

    fFunction->GetRandom2(dz, dt);

    dt *= c_light * 1.0E3; // necessary in order to make t in mm/c
    dz *= 1.0E3; // necessary in order to make z in mm

    dphi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());

    vx = 0.0;
    vy = 0.0;
    numberOfParticles = fPythia->event.size();
    for(i = 1; i < numberOfParticles; ++i)
    {
      Pythia8::Particle &particle = fPythia->event[i];

      status = particle.statusHepMC();

      if(status != 1 || !particle.isVisible() || particle.pT() <= fPTMin) continue;

      pid = particle.id();
      px = particle.px();
      py = particle.py();
      pz = particle.pz();
      e = particle.e();
      x = particle.xProd();
      y = particle.yProd();
      z = particle.zProd();
      t = particle.tProd();

      auto *candidate = factory->NewCandidate();

      candidate->PID = pid;

      candidate->Status = 1;

      pdgParticle = pdg->GetParticle(pid);
      candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge() / 3.0) : -999;
      candidate->Mass = pdgParticle ? pdgParticle->Mass() : -999.9;

      candidate->IsPU = 1;

      candidate->Momentum.SetPxPyPzE(px, py, pz, e);
      candidate->Momentum = ROOT::Math::RotationZ(dphi) * candidate->Momentum;

      x -= fInputBeamSpotX;
      y -= fInputBeamSpotY;
      candidate->Position.SetXYZT(x, y, z + dz, t + dt);
      candidate->Position = ROOT::Math::RotationZ(dphi) * candidate->Position;
      candidate->Position += ROOT::Math::XYZTVector(fOutputBeamSpotX, fOutputBeamSpotY, 0.0, 0.0);

      vx += candidate->Position.X();
      vy += candidate->Position.Y();

      fParticleOutputArray->emplace_back(*candidate);
    }

    if(numberOfParticles > 0)
    {
      vx /= numberOfParticles;
      vy /= numberOfParticles;
    }

    auto *vertex = factory->NewCandidate();
    vertex->Position.SetXYZT(vx, vy, dz, dt);
    vertex->IsPU = 1;

    fVertexOutputArray->emplace_back(*vertex);
  }
}

//------------------------------------------------------------------------------
