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

/** \class PdgCodeFilter
 *
 *  Removes particles with specific PDG codes
 *
 *  \author M. Selvaggi
 *
 */

#include "modules/PdgCodeFilter.h"

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

PdgCodeFilter::PdgCodeFilter() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

PdgCodeFilter::~PdgCodeFilter()
{
}

//------------------------------------------------------------------------------

void PdgCodeFilter::Init()
{

  ExRootConfParam param;
  Size_t i, size;

  // PT threshold
  fPTMin = GetDouble("PTMin", 0.0);

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));
  fItInputArray = fInputArray->MakeIterator();

  param = GetParam("PdgCode");
  size = param.GetSize();

  // read PdgCodes to be filtered out from the data card

  fPdgCodes.clear();
  for(i = 0; i < size; ++i)
  {
    fPdgCodes.push_back(param[i].GetInt());
  }

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "filteredParticles"));
}

//------------------------------------------------------------------------------

void PdgCodeFilter::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void PdgCodeFilter::Process()
{
  Candidate *candidate;
  Int_t pdgCode;
  Bool_t pass;
  Double_t pt;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    pdgCode = candidate->PID;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    pt = candidateMomentum.Pt();

    pass = kTRUE;

    if(pt < fPTMin) pass = kFALSE;
    if(find(fPdgCodes.begin(), fPdgCodes.end(), pdgCode) != fPdgCodes.end()) pass = kFALSE;

    if(pass) fOutputArray->Add(candidate);
  }
}

