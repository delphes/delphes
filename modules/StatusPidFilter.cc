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

/** \class StatusPidFilter
 *
 *  Removes all generated particles except electrons, muons, taus,
 *  and particles with status == 3.
 *
 *  \author J. Hirschauer - FNAL
 *
 */

#include "modules/StatusPidFilter.h"

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

namespace
{
// integer power (faster than TMath::Pow() + cast to integer)
int ipow(int base, int exp)
{
  int result = 1;
  while(exp)
  {
    if(exp & 1)
      result *= base;
    exp >>= 1;
    base *= base;
  }

  return result;
}

// standalone function to extract the i-th digit from a number (counting from 0 = rightmost, etc..)
int digit(int val, int i)
{
  int y = ipow(10, i);
  int z = val / y;
  int val2 = val / (y * 10);
  return (z - val2 * 10);
}

//  return the first two digits if this is a "fundamental" particle
//  ID = 100 is a special case (internal generator ID's are 81-100)
//  also, 101 and 102 are now used (by HepPID) for geantinos
int fundamentalID(int pdgCode)
{
  pdgCode = abs(pdgCode);
  if((digit(pdgCode, 9) == 1) && (digit(pdgCode, 8) == 0))
  {
    return 0;
  }
  if(digit(pdgCode, 2) == 0 && digit(pdgCode, 3) == 0)
  {
    return pdgCode % 10000;
  }
  else if(pdgCode <= 102)
  {
    return pdgCode;
  }
  else
  {
    return 0;
  }
}

bool hasBottom(int pdgCode)
{
  if((pdgCode / 10000000) > 0)
    return false;
  if(pdgCode <= 100)
    return false;
  if(fundamentalID(pdgCode) <= 100 && fundamentalID(pdgCode) > 0)
    return false;
  if(digit(pdgCode, 3) == 5 || digit(pdgCode, 2) == 5 || digit(pdgCode, 1) == 5)
    return true;
  return false;
}

bool isTauDaughter(int pdgCode, int M1, const TObjArray *fInputArray)
{
  //not needed, just to speed up the code - can be further refined but gives only negligible improvement:
  if(pdgCode == 15 || pdgCode < 11 || (pdgCode > 22 && pdgCode < 100) || pdgCode > 1000)
    return false;

  if(M1 < 0)
    return false;

  Candidate *mother;
  mother = static_cast<Candidate *>(fInputArray->At(M1));
  if(TMath::Abs(mother->PID) == 15)
    return true;

  return false;
}

bool isWDaughter(int M1, const TObjArray *fInputArray)
{
  if(M1 < 0) return false;

  Candidate *mother;
  mother = static_cast<Candidate *>(fInputArray->At(M1));
  if(TMath::Abs(mother->PID) == 24) return true;

  return false;
}

} // namespace

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

  // keep or remove pileup particles
  fRequireNotPileup = GetBool("RequireNotPileup", false);

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
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    status = candidate->Status;
    pdgCode = TMath::Abs(candidate->PID);

    pass = kFALSE;

    // Store all SUSY particles
    if(pdgCode >= 1000001 && pdgCode <= 1000039) pass = kTRUE;

    // hard scattering particles (first condition for Py6, second for Py8)
    if(status == 3) pass = kTRUE;
    if(status > 20 && status < 30) pass = kTRUE;

    // electrons, muons, taus and neutrinos
    if(pdgCode > 10 && pdgCode < 17) pass = kTRUE;

    // heavy quarks
    if(pdgCode == 4 || pdgCode == 5 || pdgCode == 6) pass = kTRUE;

    // Gauge bosons and other fundamental bosons
    if(pdgCode > 22 && pdgCode < 43) pass = kTRUE;

    //Stable photons
    if(pdgCode == 22 && status == 1) pass = kTRUE;

    // logic ported from HepPDF: http://lcgapp.cern.ch/project/simu/HepPDT/HepPDT.2.05.02/html/ParticleID_8cc-source.html#l00081
    bool is_b_hadron = hasBottom(pdgCode);
    bool is_b_quark = (pdgCode == 5);

    bool is_tau_daughter = isTauDaughter(pdgCode, candidate->M1, fInputArray);

    if(is_b_hadron)
      pass = kTRUE;

    if(is_tau_daughter)
      pass = kTRUE;

    bool is_W_daughter = isWDaughter(candidate->M1, fInputArray);
    if(is_W_daughter)
      pass = kTRUE;

    // fPTMin not applied to b_hadrons / b_quarks to allow for b-enriched sample stitching
    // fPTMin not applied to tau decay products to allow visible-tau four momentum determination
    if(!pass || (candidate->Momentum.Pt() < fPTMin && !(is_b_hadron || is_b_quark || is_tau_daughter || is_W_daughter))) continue;

    // not pileup particles
    if(fRequireNotPileup && (candidate->IsPU > 0)) continue;

    fOutputArray->Add(candidate);
  }
}
