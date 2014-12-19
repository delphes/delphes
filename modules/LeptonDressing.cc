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


/** \class LeptonDressing
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
