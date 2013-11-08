/** \class PileUpMerger
 *
 *  Merges particles from pile-up sample into event
 *
 *
 *  $Date: 2013-02-12 15:13:59 +0100 (Tue, 12 Feb 2013) $
 *  $Revision: 907 $
 *
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "modules/PileUpMerger.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
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
  fReader(0), fItInputArray(0)
{
}

//------------------------------------------------------------------------------

PileUpMerger::~PileUpMerger()
{
}

//------------------------------------------------------------------------------

void PileUpMerger::Init()
{
  const char *fileName;

  fPileUpDistribution = GetInt("PileUpDistribution", 0);

  fMeanPileUp  = GetDouble("MeanPileUp", 10);
  fZVertexSpread = GetDouble("ZVertexSpread", 0.05)*1.0E3;

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
  Double_t dz, dphi;
  Int_t numberOfEvents, event;
  Long64_t allEntries, entry;
  Candidate *candidate;
  DelphesFactory *factory;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    fParticleOutputArray->Add(candidate);
  }

  factory = GetFactory();

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

    dz = gRandom->Gaus(0.0, fZVertexSpread);
    dphi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());

    candidate = factory->NewCandidate();
    candidate->Position.SetXYZT(0.0, 0.0, dz, 0.0);
    fVertexOutputArray->Add(candidate);

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

      candidate->Position.SetXYZT(x, y, z + dz, t);
      candidate->Position.RotateZ(dphi);

      fParticleOutputArray->Add(candidate);
    }
  }
}

//------------------------------------------------------------------------------

