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

#include "modules/PileUpMerger.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesTF2.h"
#include "classes/DelphesPileUpReader.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

PileUpMerger::PileUpMerger() :
  fFunction(0), fReader(0), fItInputArray(0)
{
  fFunction = new DelphesTF2;
}


//------------------------------------------------------------------------------

PileUpMerger::~PileUpMerger()
{
  delete fFunction;
}

//------------------------------------------------------------------------------

void PileUpMerger::Init()
{
  const char *fileName;

  fPileUpDistribution = GetInt("PileUpDistribution", 0);

  fMeanPileUp  = GetDouble("MeanPileUp", 10);

  fZVertexSpread = GetDouble("ZVertexSpread", 0.15);
  fTVertexSpread = GetDouble("TVertexSpread", 1.5E-09);

  // read vertex smearing formula

  fFunction->Compile(GetString("VertexDistributionFormula", "0.0"));
  fFunction->SetRange(-fZVertexSpread, -fTVertexSpread, fZVertexSpread, fTVertexSpread);

  fileName = GetString("PileUpFile", "MinBias.pileup");
  fReader = new DelphesPileUpReader(fileName);

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output arrays
  fParticleOutputArray = ExportArray(GetString("ParticleOutputArray", "stableParticles"));
  fVertexOutputArray = ExportArray(GetString("VertexOutputArray", "vertices"));
}

//------------------------------------------------------------------------------

void PileUpMerger::Finish()
{
  if(fReader) delete fReader;
}

//------------------------------------------------------------------------------

void PileUpMerger::Process()
{
  TDatabasePDG *pdg = TDatabasePDG::Instance();
  TParticlePDG *pdgParticle;
  Int_t pid;
  Float_t x, y, z, t;
  Float_t px, py, pz, e;
  Double_t dz, dphi, dt;
  Int_t numberOfEvents, event;
  Long64_t allEntries, entry;
  Candidate *candidate, *vertexcandidate;
  DelphesFactory *factory;

  const Double_t c_light = 2.99792458E8;

  fItInputArray->Reset();

  // --- Deal with Primary vertex first  ------

  fFunction->GetRandom2(dz, dt);

  dt *= c_light*1.0E3; // necessary in order to make t in mm/c
  dz *= 1.0E3; // necessary in order to make z in mm

  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    z = candidate->Position.Z();
    t = candidate->Position.T();
    candidate->Position.SetZ(z + dz);
    candidate->Position.SetT(t + dt);
    fParticleOutputArray->Add(candidate);
  }

  factory = GetFactory();

  vertexcandidate = factory->NewCandidate();
  vertexcandidate->Position.SetXYZT(0.0, 0.0, dz, dt);
  fVertexOutputArray->Add(vertexcandidate);

  // --- Then with pile-up vertices  ------

  switch(fPileUpDistribution)
  {
    case 0:
      numberOfEvents = gRandom->Poisson(fMeanPileUp);
      break;
    case 1:
      numberOfEvents = gRandom->Integer(2*fMeanPileUp + 1);
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
      entry = TMath::Nint(gRandom->Rndm()*allEntries);
    }
    while(entry >= allEntries);

    fReader->ReadEntry(entry);

   // --- Pile-up vertex smearing

    fFunction->GetRandom2(dz, dt);

    dt *= c_light*1.0E3; // necessary in order to make t in mm/c
    dz *= 1.0E3; // necessary in order to make z in mm

    dphi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());

    vertexcandidate = factory->NewCandidate();
    vertexcandidate->Position.SetXYZT(0.0, 0.0, dz, dt);
    vertexcandidate->IsPU = 1;

    fVertexOutputArray->Add(vertexcandidate);

    while(fReader->ReadParticle(pid, x, y, z, t, px, py, pz, e))
    {
      candidate = factory->NewCandidate();

      candidate->PID = pid;

      candidate->Status = 1;

      pdgParticle = pdg->GetParticle(pid);
      candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge()/3.0) : -999;
      candidate->Mass = pdgParticle ? pdgParticle->Mass() : -999.9;

      candidate->IsPU = 1;

      candidate->Momentum.SetPxPyPzE(px, py, pz, e);
      candidate->Momentum.RotateZ(dphi);

      candidate->Position.SetXYZT(x, y, z + dz, t + dt);
      candidate->Position.RotateZ(dphi);

      fParticleOutputArray->Add(candidate);
    }
  }
}

//------------------------------------------------------------------------------

