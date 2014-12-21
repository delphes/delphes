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


/** \class Merger
 *
 *  Merges multiple input arrays into one output array
 *  and sums transverse momenta of all input objects.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/Merger.h"

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

Merger::Merger()
{
}

//------------------------------------------------------------------------------

Merger::~Merger()
{
}

//------------------------------------------------------------------------------

void Merger::Init()
{
  // import arrays with output from other modules

  ExRootConfParam param = GetParam("InputArray");
  Long_t i, size;
  const TObjArray *array;
  TIterator *iterator;

  size = param.GetSize();
  for(i = 0; i < size; ++i)
  {
    array = ImportArray(param[i].GetString());
    iterator = array->MakeIterator();

    fInputList.push_back(iterator);
  }

  // create output arrays

  fOutputArray = ExportArray(GetString("OutputArray", "candidates"));

  fMomentumOutputArray = ExportArray(GetString("MomentumOutputArray", "momentum"));
  
  fEnergyOutputArray = ExportArray(GetString("EnergyOutputArray", "energy"));
}

//------------------------------------------------------------------------------

void Merger::Finish()
{
  vector< TIterator * >::iterator itInputList;
  TIterator *iterator;

  for(itInputList = fInputList.begin(); itInputList != fInputList.end(); ++itInputList)
  {
    iterator = *itInputList;
    if(iterator) delete iterator;
  }
}

//------------------------------------------------------------------------------

void Merger::Process()
{
  Candidate *candidate;
  TLorentzVector momentum;
  Double_t sumPT, sumE;  
  vector< TIterator * >::iterator itInputList;
  TIterator *iterator;

  DelphesFactory *factory = GetFactory();
  
  momentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
  sumPT = 0;
  sumE = 0;

  // loop over all input arrays
  for(itInputList = fInputList.begin(); itInputList != fInputList.end(); ++itInputList)
  {
    iterator = *itInputList;

    // loop over all candidates
    iterator->Reset();
    while((candidate = static_cast<Candidate*>(iterator->Next())))
    {
      const TLorentzVector &candidateMomentum = candidate->Momentum;

      momentum += candidateMomentum;
      sumPT += candidateMomentum.Pt();
      sumE += candidateMomentum.E();

      fOutputArray->Add(candidate);
    }
  }

  candidate = factory->NewCandidate();
  
  candidate->Position.SetXYZT(0.0, 0.0, 0.0, 0.0);
  candidate->Momentum = momentum;
  
  fMomentumOutputArray->Add(candidate);

  candidate = factory->NewCandidate();
  
  candidate->Position.SetXYZT(0.0, 0.0, 0.0, 0.0);
  candidate->Momentum.SetPtEtaPhiE(sumPT, 0.0, 0.0, sumE);
  
  fEnergyOutputArray->Add(candidate);
}

//------------------------------------------------------------------------------
