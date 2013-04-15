
/** \class LeptonDressing
 *
 *
 *  
 *
 *
 *  \author P. Demin && A. Mertens - UCL, Louvain-la-Neuve
 *
 */

#include "modules/LeptonDressing.h"

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

LeptonDressing::LeptonDressing() :
 fItDressingInputArray(0), fItCandidateInputArray(0)
{
}

//------------------------------------------------------------------------------

LeptonDressing::~LeptonDressing()
{
}

//------------------------------------------------------------------------------

void LeptonDressing::Init()
{
  fDeltaR = GetDouble("DeltaRMax", 0.4);

  // import input array(s)

  fDressingInputArray = ImportArray(GetString("DressingInputArray", "Calorimeter/photons"));
  fItDressingInputArray = fDressingInputArray->MakeIterator();
  
  fCandidateInputArray = ImportArray(GetString("CandidateInputArray", "UniqueObjectFinder/electrons"));
  fItCandidateInputArray = fCandidateInputArray->MakeIterator();
  
  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "electrons"));
}

//------------------------------------------------------------------------------

void LeptonDressing::Finish()
{
  if(fItCandidateInputArray) delete fItCandidateInputArray;
  if(fItDressingInputArray) delete fItDressingInputArray;
}

//------------------------------------------------------------------------------

void LeptonDressing::Process()
{
  Candidate *candidate, *dressing, *mother;
  TLorentzVector momentum;
  
  // loop over all input candidate
  fItCandidateInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItCandidateInputArray->Next())))
  {
    const TLorentzVector &candidateMomentum = candidate->Momentum;

    // loop over all input tracks
    fItDressingInputArray->Reset();
    momentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
    while((dressing = static_cast<Candidate*>(fItDressingInputArray->Next())))
    {
      const TLorentzVector &dressingMomentum = dressing->Momentum;
      if (dressingMomentum.Pt() > 0.1)
      {
        if(candidateMomentum.DeltaR(dressingMomentum) <= fDeltaR)
        {
          momentum += dressingMomentum;
        }
      }
    }

    mother = candidate;
    candidate = static_cast<Candidate*>(candidate->Clone());

    candidate->Momentum += momentum;
    candidate->AddCandidate(mother);
    
    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
