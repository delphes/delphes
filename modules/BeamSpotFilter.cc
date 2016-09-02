/** \class BeamSpotFilter
 *
 *  Extracts beam spot 
 *
 *  \author Michele Selvaggi
 *
 */


#include "modules/BeamSpotFilter.h"

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

BeamSpotFilter::BeamSpotFilter() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

BeamSpotFilter::~BeamSpotFilter()
{
}

//------------------------------------------------------------------------------

void BeamSpotFilter::Init()
{
 
  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "filteredParticles"));
}

//------------------------------------------------------------------------------

void BeamSpotFilter::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void BeamSpotFilter::Process()
{
  Candidate *candidate;
  Bool_t passed = false;
  
  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())) && !passed)
  {
    if(candidate->IsPU == 0) passed = true;
    fOutputArray->Add(candidate);
  }
  
}

 
