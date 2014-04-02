
/** \class Hector
 *
 *  Propagates candidates using Hector library.
 *
 *  $Date$
 *  $Revision$
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/Hector.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

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

#include "Hector/H_BeamLine.h"
#include "Hector/H_RecRPObject.h"
#include "Hector/H_BeamParticle.h"

using namespace std;

//------------------------------------------------------------------------------

Hector::Hector() :
  fBeamLine(0), fItInputArray(0)
{
}

//------------------------------------------------------------------------------

Hector::~Hector()
{
}

//------------------------------------------------------------------------------

void Hector::Init()
{
  // read Hector parameters

  fDirection = GetInt("Direction", 1);
  fBeamLineLength = GetDouble("BeamLineLength", 430.0);
  fDistance = GetDouble("Distance", 420.0);
  fSigmaE = GetDouble("SigmaE", 0.0);
  fSigmaX = GetDouble("SigmaX", 0.0);
  fSigmaY = GetDouble("SigmaY", 0.0);
  fEtaMin = GetDouble("EtaMin", 5.0);

  fBeamLine = new H_BeamLine(fDirection, fBeamLineLength + 0.1);
  fBeamLine->fill(GetString("BeamLineFile", "examples/LHCB1IR5_5TeV.tfs"), fDirection, "IP5" );
  fBeamLine->offsetElements(120, -0.097*fDirection);
  fBeamLine->calcMatrix();

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void Hector::Finish()
{
  if(fItInputArray) delete fItInputArray;
  if(fBeamLine) delete fBeamLine;
}

//------------------------------------------------------------------------------

void Hector::Process()
{
  Candidate *candidate, *mother;
  Double_t pz;
  Double_t x, y, z, tx, ty;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    pz = candidateMomentum.Pz();

    if(candidateMomentum.Eta() <= fEtaMin || fDirection*pz <= 0.0) continue;

    x = 1.0E3 * candidatePosition.X();
    y = 1.0E3 * candidatePosition.Y();
    z = 1.0E-2 * candidatePosition.Z();

//    tx = 1.0E6 * TMath::ATan(candidateMomentum.Px()/pz);
//    ty = 1.0E6 * TMath::ATan(candidateMomentum.Py()/pz);

    tx = 0.0;
    ty = 0.0;

    H_BeamParticle particle(candidate->Mass, candidate->Charge);
    particle.set4Momentum(candidateMomentum);
    particle.setPosition(x, y, tx, ty, z);

    particle.smearAng(fSigmaX, fSigmaY, gRandom);
    particle.smearE(fSigmaE, gRandom);

    particle.computePath(fBeamLine);

    if(particle.stopped(fBeamLine)) continue;

    particle.propagate(fDistance);

    mother = candidate;
    candidate = static_cast<Candidate*>(candidate->Clone());
    candidate->Position.SetXYZT(particle.getX(), particle.getY(), fDistance, 0.0);
    candidate->Momentum.SetPxPyPzE(particle.getTX(), particle.getTY(), 0.0, particle.getE());
    candidate->AddCandidate(mother);

    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
