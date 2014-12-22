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
#include "classes/DelphesFormula.h"
#include "classes/DelphesPileUpReader.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "Pythia.h"

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

PileUpMergerPythia8::PileUpMergerPythia8() :
  fPythia(0), fItInputArray(0)
{
}

//------------------------------------------------------------------------------

PileUpMergerPythia8::~PileUpMergerPythia8()
{
}

//------------------------------------------------------------------------------

void PileUpMergerPythia8::Init()
{
  const char *fileName;

  fMeanPileUp  = GetDouble("MeanPileUp", 10);
  fZVertexSpread = GetDouble("ZVertexSpread", 0.05)*1.0E3;

  fPTMin = GetDouble("PTMin", 0.0);

  fileName = GetString("ConfigFile", "MinBias.cmnd");
  fPythia = new Pythia8::Pythia();
  fPythia->readFile(fileName);

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output arrays
  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
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
  Float_t x, y, z, t;
  Float_t px, py, pz, e;
  Double_t dz, dphi;
  Int_t poisson, event, i;
  Candidate *candidate;
  DelphesFactory *factory;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    fOutputArray->Add(candidate);
  }

  factory = GetFactory();

  poisson = gRandom->Poisson(fMeanPileUp);

  for(event = 0; event < poisson; ++event)
  {
    while(!fPythia->next());

    dz = gRandom->Gaus(0.0, fZVertexSpread);
    dphi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());

    for(i = 1; i < fPythia->event.size(); ++i)
    {
      Pythia8::Particle &particle = fPythia->event[i];

      status = particle.statusHepMC();

      if(status != 1 || !particle.isVisible() || particle.pT() <= fPTMin) continue;

      pid = particle.id();
      px = particle.px(); py = particle.py(); pz = particle.pz(); e = particle.e();
      x = particle.xProd(); y = particle.yProd(); z = particle.zProd(); t = particle.tProd();

      candidate = factory->NewCandidate();

      candidate->PID = pid;

      candidate->Status = status;

      pdgParticle = pdg->GetParticle(pid);
      candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge()/3.0) : -999;
      candidate->Mass = pdgParticle ? pdgParticle->Mass() : -999.9;

      candidate->IsPU = 1;

      candidate->Momentum.SetPxPyPzE(px, py, pz, e);
      candidate->Momentum.RotateZ(dphi);

      candidate->Position.SetXYZT(x, y, z + dz, t);
      candidate->Position.RotateZ(dphi);

      fOutputArray->Add(candidate);
    }
  }
}

//------------------------------------------------------------------------------

