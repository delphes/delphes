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


/** \class PileUpJetID
 *
 *  CMS PileUp Jet ID Variables, based on http://cds.cern.ch/record/1581583
 *
 *  \author S. Zenz, December 2013
 *
 */

#include "modules/PileUpJetID.h"

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

PileUpJetID::PileUpJetID() :
  fItJetInputArray(0),fTrackInputArray(0),fNeutralInputArray(0),fItVertexInputArray(0)
{

}

//------------------------------------------------------------------------------

PileUpJetID::~PileUpJetID()
{

}

//------------------------------------------------------------------------------

void PileUpJetID::Init()
{
  fJetPTMin = GetDouble("JetPTMin", 20.0);
  fParameterR = GetDouble("ParameterR", 0.5);
  fUseConstituents = GetInt("UseConstituents", 0);

  fAverageEachTower = false; // for timing

  // import input array(s)

  fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();

  fTrackInputArray = ImportArray(GetString("TrackInputArray", "Calorimeter/eflowTracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();

  fNeutralInputArray = ImportArray(GetString("NeutralInputArray", "Calorimeter/eflowTowers"));
  fItNeutralInputArray = fNeutralInputArray->MakeIterator();

  fVertexInputArray = ImportArray(GetString("VertexInputArray", "PileUpMerger/vertices"));
  fItVertexInputArray = fVertexInputArray->MakeIterator();

  fZVertexResolution  = GetDouble("ZVertexResolution", 0.005)*1.0E3;

  // create output array(s)

  fOutputArray = ExportArray(GetString("OutputArray", "jets"));
}

//------------------------------------------------------------------------------

void PileUpJetID::Finish()
{
  if(fItJetInputArray) delete fItJetInputArray;
  if(fItTrackInputArray) delete fItTrackInputArray;
  if(fItNeutralInputArray) delete fItNeutralInputArray;
  if(fItVertexInputArray) delete fItVertexInputArray;
}

//------------------------------------------------------------------------------

void PileUpJetID::Process()
{
  Candidate *candidate, *constituent;
  TLorentzVector momentum, area;
  Int_t i, nc, nn;
  Double_t sumpt, sumptch, sumptchpv, sumptchpu, sumdrsqptsq, sumptsq;
  Double_t dr, pt, pt_ann[5];
  Double_t zvtx = 0.0;

  Candidate *track;

  // find z position of primary vertex

  fItVertexInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItVertexInputArray->Next())))
  {
    if(!candidate->IsPU)
    {
      zvtx = candidate->Position.Z();
      break;
    }
  }

  // loop over all input candidates
  fItJetInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItJetInputArray->Next())))
  {
    momentum = candidate->Momentum;
    area = candidate->Area;

    sumpt = 0.0;
    sumptch = 0.0;
    sumptchpv = 0.0;
    sumptchpu = 0.0;
    sumdrsqptsq = 0.0;
    sumptsq = 0.0;
    nc = 0;
    nn = 0;

    for(i = 0; i < 5; ++i)
    {
      pt_ann[i] = 0.0;
    }

    if(fUseConstituents)
    {
      TIter itConstituents(candidate->GetCandidates());
      while((constituent = static_cast<Candidate*>(itConstituents.Next())))
      {
        pt = constituent->Momentum.Pt();
        dr = candidate->Momentum.DeltaR(constituent->Momentum);
        sumpt += pt;
        sumdrsqptsq += dr*dr*pt*pt;
        sumptsq += pt*pt;
        if(constituent->Charge == 0)
        {
          // neutrals
          ++nn;
        }
        else
        {
          // charged
          if(constituent->IsPU && TMath::Abs(constituent->Position.Z()-zvtx) > fZVertexResolution)
          {
            sumptchpu += pt;
          }
          else
          {
            sumptchpv += pt;
          }
          sumptch += pt;
          ++nc;
        }
        for(i = 0; i < 5; ++i)
        {
          if(dr > 0.1*i && dr < 0.1*(i + 1))
          {
            pt_ann[i] += pt;
          }
        }
      }
    }
    else
    {
      // Not using constituents, using dr
      fItTrackInputArray->Reset();
      while((track = static_cast<Candidate*>(fItTrackInputArray->Next())))
      {
        if(track->Momentum.DeltaR(candidate->Momentum) < fParameterR)
        {
          pt = track->Momentum.Pt();
          sumpt += pt;
          sumptch += pt;
          if(track->IsPU && TMath::Abs(track->Position.Z()-zvtx) > fZVertexResolution)
          {
            sumptchpu += pt;
          }
          else
          {
            sumptchpv += pt;
          }
          dr = candidate->Momentum.DeltaR(track->Momentum);
          sumdrsqptsq += dr*dr*pt*pt;
          sumptsq += pt*pt;
          nc++;
          for(i = 0; i < 5; ++i)
          {
            if(dr > 0.1*i && dr < 0.1*(i + 1))
            {
              pt_ann[i] += pt;
            }
          }
        }
      }

      fItNeutralInputArray->Reset();
      while ((constituent = static_cast<Candidate*>(fItNeutralInputArray->Next())))
      {
        if(constituent->Momentum.DeltaR(candidate->Momentum) < fParameterR)
        {
          pt = constituent->Momentum.Pt();
          sumpt += pt;
          dr = candidate->Momentum.DeltaR(constituent->Momentum);
          sumdrsqptsq += dr*dr*pt*pt;
          sumptsq += pt*pt;
          nn++;
          for(i = 0; i < 5; ++i)
          {
            if(dr > 0.1*i && dr < 0.1*(i + 1))
            {
              pt_ann[i] += pt;
            }
          }
        }
      }
    }

    if(sumptch > 0.0)
    {
      candidate->Beta = sumptchpu/sumptch;
      candidate->BetaStar = sumptchpv/sumptch;
    }
    else
    {
      candidate->Beta = -999.0;
      candidate->BetaStar = -999.0;
    }
    if(sumptsq > 0.0)
    {
      candidate->MeanSqDeltaR = sumdrsqptsq/sumptsq;
    }
    else
    {
      candidate->MeanSqDeltaR = -999.0;
    }
    candidate->NCharged = nc;
    candidate->NNeutrals = nn;
    if(sumpt > 0.0)
    {
      candidate->PTD = TMath::Sqrt(sumptsq) / sumpt;
      for(i = 0; i < 5; ++i)
      {
        candidate->FracPt[i] = pt_ann[i]/sumpt;
      }
    }
    else
    {
      candidate->PTD = -999.0;
      for(i = 0; i < 5; ++i)
      {
        candidate->FracPt[i] = -999.0;
      }
    }

    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------

