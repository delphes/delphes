
/** \class AngularSmearing
 *
 *  Performs transverse angular resolution smearing.
 *
 *  $Date: 2014-06-17 16:58:53 +0100  $
 *  
 *
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "modules/AngularSmearing.h"

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

AngularSmearing::AngularSmearing() :
  fFormulaEta(0), fFormulaPhi(0), fItInputArray(0)
{
  fFormulaEta = new DelphesFormula;
  fFormulaPhi = new DelphesFormula;
}

//------------------------------------------------------------------------------

AngularSmearing::~AngularSmearing()
{
  if(fFormulaEta) delete fFormulaEta;
  if(fFormulaPhi) delete fFormulaPhi;
}

//------------------------------------------------------------------------------

void AngularSmearing::Init()
{
  // read resolution formula

  fFormulaEta->Compile(GetString("EtaResolutionFormula", "0.0"));
  fFormulaPhi->Compile(GetString("PhiResolutionFormula", "0.0"));


  // import input array

  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void AngularSmearing::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void AngularSmearing::Process()
{
  Candidate *candidate, *mother;
  Double_t pt, eta, phi;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    eta = candidatePosition.Eta();
    phi = candidatePosition.Phi();
    pt = candidateMomentum.Pt();

    // apply smearing formula for eta,phi

    eta = gRandom->Gaus(eta, fFormulaEta->Eval(pt, eta));
    phi = gRandom->Gaus(phi, fFormulaPhi->Eval(pt, eta));
    
    if(pt <= 0.0) continue;

    mother = candidate;
    candidate = static_cast<Candidate*>(candidate->Clone());
    eta = candidateMomentum.Eta();
    phi = candidateMomentum.Phi();
    candidate->Momentum.SetPtEtaPhiE(pt, eta, phi, pt*TMath::CosH(eta));
    candidate->AddCandidate(mother);
        
    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
