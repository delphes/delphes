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

/** \class TrackPileUpSubtractor
 *
 *  Subtract pile-up contribution from tracks.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *  \author J. T. Offermann - Brown University
 * (J.T.O. adding the sumpt2 vertex option)
 *
 */

#include "modules/TrackPileUpSubtractor.h"

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

TrackPileUpSubtractor::TrackPileUpSubtractor() :
  fFormula(0)
{
  fFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

TrackPileUpSubtractor::~TrackPileUpSubtractor()
{
  if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void TrackPileUpSubtractor::Init()
{
  // import input array

  fVertexInputArray = ImportArray(GetString("VertexInputArray", "PileUpMerger/vertices"));
  fItVertexInputArray = fVertexInputArray->MakeIterator();

  // read resolution formula in m
  fFormula->Compile(GetString("ZVertexResolution", "0.001"));

  fPTMin = GetDouble("PTMin", 0.);

  // determine mode -- do we define the primary vertex using truth-level info (Candidate::IsPU),
  // or do we define it as the vertex with the largest sum(pt^2) of charged particles?
  fMode = GetInt("PrimaryVertexMode",0);

  // import arrays with output from other modules

  ExRootConfParam param = GetParam("InputArray");
  Long_t i, size;
  const TObjArray *array;
  TIterator *iterator;

  size = param.GetSize();
  for(i = 0; i < size / 2; ++i)
  {
    array = ImportArray(param[i * 2].GetString());
    iterator = array->MakeIterator();

    fInputMap[iterator] = ExportArray(param[i * 2 + 1].GetString());
  }
}

//------------------------------------------------------------------------------

void TrackPileUpSubtractor::Finish()
{
  map<TIterator *, TObjArray *>::iterator itInputMap;
  TIterator *iterator;

  for(itInputMap = fInputMap.begin(); itInputMap != fInputMap.end(); ++itInputMap)
  {
    iterator = itInputMap->first;

    if(iterator) delete iterator;
  }

  if(fItVertexInputArray) delete fItVertexInputArray;
}

//------------------------------------------------------------------------------

void TrackPileUpSubtractor::Process()
{
  Candidate *candidate, *particle;
  map<TIterator *, TObjArray *>::iterator itInputMap;
  TIterator *iterator;
  TObjArray *array;
  Double_t z, zvtx = 0;
  Double_t pt, eta, phi, e;
  Double_t sumpt2, sumpt2_max = -999.;

  // find z position of primary vertex

  fItVertexInputArray->Reset();

  if(fMode == 1)
  {
    while((candidate = static_cast<Candidate *>(fItVertexInputArray->Next())))
    {
      sumpt2 = candidate->SumPT2;
      if(sumpt2 == 0.) sumpt2 = candidate->GenSumPT2; // in some cases, only one or the other is filled (they default to zero)
      if(sumpt2 > sumpt2_max)
      {
        sumpt2_max = sumpt2;
        zvtx = candidate->Position.Z();
      }
    }
  }
  else // for now, only 2 options for fMode, hence if/else should probably be OK.
  {
    while((candidate = static_cast<Candidate *>(fItVertexInputArray->Next())))
    {
      if(!candidate->IsPU)
      {
        zvtx = candidate->Position.Z();
        // break;
      }
    }
}


  // loop over all input arrays
  for(itInputMap = fInputMap.begin(); itInputMap != fInputMap.end(); ++itInputMap)
  {
    iterator = itInputMap->first;
    array = itInputMap->second;

    // loop over all candidates
    iterator->Reset();
    while((candidate = static_cast<Candidate *>(iterator->Next())))
    {
      particle = static_cast<Candidate *>(candidate->GetCandidates()->At(0));
      const TLorentzVector &candidateMomentum = particle->Momentum;

      eta = candidateMomentum.Eta();
      pt = candidateMomentum.Pt();
      phi = candidateMomentum.Phi();
      e = candidateMomentum.E();

      z = particle->Position.Z();

      // apply pile-up subtraction
      // assume perfect pile-up subtraction for tracks outside fZVertexResolution

      if(candidate->Charge != 0 && candidate->IsPU && TMath::Abs(z - zvtx) > fFormula->Eval(pt, eta, phi, e) * 1.0e3)
      {
        candidate->IsRecoPU = 1;
      }
      else
      {
        candidate->IsRecoPU = 0;
        if(candidate->Momentum.Pt() > fPTMin) array->Add(candidate);
      }
    }
  }
}
