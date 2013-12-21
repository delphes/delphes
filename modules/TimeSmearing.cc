
/** \class TimeSmearing
 *
 *  Performs transverse momentum resolution smearing.
 *
 *  $Date: 2013-02-13 16:58:53 +0100 (Wed, 13 Feb 2013) $
 *  $Revision: 911 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/TimeSmearing.h"

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

TimeSmearing::TimeSmearing() :
fItInputArray(0)
{
}

//------------------------------------------------------------------------------

TimeSmearing::~TimeSmearing()
{
}

//------------------------------------------------------------------------------

void TimeSmearing::Init()
{
  // read resolution formula

  fTimeResolution = GetDouble("TimeResolution", 1.0E-10);
  // import input array

  fInputArray = ImportArray(GetString("InputArray", "MuonMomentumSmearing/muons"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "muons"));
}

//------------------------------------------------------------------------------

void TimeSmearing::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void TimeSmearing::Process()
{
  Candidate *candidate, *mother;
  Double_t t;
  const Double_t c_light = 2.99792458E8;
  
  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    t = candidatePosition.T()*1.0E-3/c_light;
    
    // apply smearing formula
    t = gRandom->Gaus(t, fTimeResolution);
   
    mother = candidate;
    candidate = static_cast<Candidate*>(candidate->Clone());
    candidate->Position.SetT(t*1.0E3*c_light);
    
    candidate->AddCandidate(mother);
        
    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
