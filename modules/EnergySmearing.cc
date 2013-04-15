
/** \class EnergySmearing
 *
 *  Performs energy resolution smearing.
 *
 *  $Date$
 *  $Revision$
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/EnergySmearing.h"

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

using namespace std;

//------------------------------------------------------------------------------

EnergySmearing::EnergySmearing() :
  fFormula(0), fItInputArray(0)
{
  fFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

EnergySmearing::~EnergySmearing()
{
  if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void EnergySmearing::Init()
{
  // read resolution formula

  fFormula->Compile(GetString("ResolutionFormula", "0.0"));

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void EnergySmearing::Finish()
{  
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void EnergySmearing::Process()
{
  Candidate *candidate, *mother;
  Double_t energy, eta, phi;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    eta = candidatePosition.Eta();
    phi = candidatePosition.Phi();
    energy = candidateMomentum.E();
 
    // apply smearing formula
    energy = gRandom->Gaus(energy, fFormula->Eval(0.0, eta, 0.0, energy));
     
    if(energy <= 0.0) continue;
 
    mother = candidate;
    candidate = static_cast<Candidate*>(candidate->Clone());
    eta = candidateMomentum.Eta();
    phi = candidateMomentum.Phi();
    candidate->Momentum.SetPtEtaPhiE(energy/TMath::CosH(eta), eta, phi, energy);
    candidate->AddCandidate(mother);
 
    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
