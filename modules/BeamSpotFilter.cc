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

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

BeamSpotFilter::BeamSpotFilter() {}

//------------------------------------------------------------------------------

BeamSpotFilter::~BeamSpotFilter() {}

//------------------------------------------------------------------------------

void BeamSpotFilter::Init()
{

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));
  fItInputArray.reset(fInputArray->MakeIterator());

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "filteredParticles"));
}

//------------------------------------------------------------------------------

void BeamSpotFilter::Finish() {}

//------------------------------------------------------------------------------

void BeamSpotFilter::Process()
{
  Candidate *candidate = nullptr;
  Bool_t passed = false;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())) && !passed)
  {
    if(candidate->IsPU == 0) passed = true;
    fOutputArray->Add(candidate);
  }
}
