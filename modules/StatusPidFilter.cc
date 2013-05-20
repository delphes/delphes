/** \class StatusPidFilter
 *
 *  Removes all generated particles except electrons, muons, taus,
 *  and particles with status == 3.
 *
 *  $Date$
 *  $Revision$
 *
 *
 *  \author J. Hirschauer - FNAL
 *
 */

#include "modules/StatusPidFilter.h"

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

StatusPidFilter::StatusPidFilter() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

StatusPidFilter::~StatusPidFilter()
{
}

//------------------------------------------------------------------------------

void StatusPidFilter::Init()
{
  // PT threshold

  fPTMin = GetDouble("PTMin", 0.5);

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "filteredParticles"));
}

//------------------------------------------------------------------------------

void StatusPidFilter::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void StatusPidFilter::Process()
{
  Candidate *candidate;
  Int_t status, pdgCode;
  Bool_t pass;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    status = candidate->Status;
    pdgCode = TMath::Abs(candidate->PID);

    pass = kFALSE;

    // status == 3
    if(status == 3) pass = kTRUE;

    // electrons, muons, taus and neutrinos
    if(pdgCode > 10 && pdgCode < 17) pass = kTRUE;

    // heavy quarks
    if(pdgCode == 5 || pdgCode == 6) pass = kTRUE;

    // Gauge bosons and other fundamental bosons
    if(pdgCode > 22 && pdgCode < 43) pass = kTRUE;

    if(!pass || candidate->Momentum.Pt() <= fPTMin) continue;

    fOutputArray->Add(candidate);
  }
}

