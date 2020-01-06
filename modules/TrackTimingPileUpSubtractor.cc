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
 *
 */

#include "modules/TrackTimingPileUpSubtractor.h"

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

TrackTimingPileUpSubtractor::TrackTimingPileUpSubtractor() :
  fFormula(0)
{
  fFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

TrackTimingPileUpSubtractor::~TrackTimingPileUpSubtractor()
{
  if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void TrackTimingPileUpSubtractor::Init()
{
  // import input array

  fVertexInputArray = ImportArray(GetString("VertexInputArray", "PileUpMerger/vertices"));
  fItVertexInputArray = fVertexInputArray->MakeIterator();

  // read resolution formula in m
  fZVertexResolution = GetDouble("ZVertexResolution", 3);
  fTVertexResolution = GetDouble("TVertexResolution", 3);

  fPTMin = GetDouble("PTMin", 0.);

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

void TrackTimingPileUpSubtractor::Finish()
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

void TrackTimingPileUpSubtractor::Process()
{
  Candidate *candidate, *particle;
  map<TIterator *, TObjArray *>::iterator itInputMap;
  TIterator *iterator;
  TObjArray *array;
  Double_t z, zvtx = 0;
  Double_t z_err, zvtx_err = 0;
  Double_t t, tvtx = 0;
  Double_t t_err, tvtx_err = 0;
  Double_t sumPTSquare = 0;
  Double_t tempPTSquare = 0;
  Double_t pt, eta, phi, e;
  Double_t distanceCharged, distanceNeutral = 0;

  // find z position of primary vertex

  fItVertexInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItVertexInputArray->Next())))
  {
    tempPTSquare = candidate->SumPT2;
    if(tempPTSquare > sumPTSquare)
    {
      sumPTSquare = tempPTSquare;
      zvtx = candidate->Position.Z();
      zvtx_err = candidate->PositionError.Z();
      tvtx = candidate->Position.T();
      tvtx_err = candidate->PositionError.T();
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
      z_err = particle->PositionError.Z();
      t = particle->InitialPosition.T();
      t_err = particle->PositionError.T();

      distanceCharged = pow((zvtx - z),2)/pow((zvtx_err - z_err),2) + pow((tvtx - t),2)/pow((tvtx_err - t_err),2);
      distanceNeutral = pow((tvtx - t),2)/pow((tvtx_err - t_err),2);

      if(candidate->Charge != 0 && distanceCharged < fZVertexResolution)
      {
        candidate->IsRecoPU = 1;
      }
      else if(candidate->Charge == 0 && distanceNeutral < fTVertexResolution)
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
