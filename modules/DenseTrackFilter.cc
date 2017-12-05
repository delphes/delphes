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


/** \class DenseTrackFilter
 *
 *  Applies reconstruction inefficiency on tracks in dense environment
 *
 *  \author M. Selvaggi - CERN
 *
 */

#include "modules/DenseTrackFilter.h"

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

DenseTrackFilter::DenseTrackFilter() :
  fItTrackInputArray(0),   fItDenseChargedInputArray(0)
{
}

//------------------------------------------------------------------------------

DenseTrackFilter::~DenseTrackFilter()
{
}

//------------------------------------------------------------------------------

void DenseTrackFilter::Init()
{
  ExRootConfParam param, paramEtaBins, paramPhiBins, paramFractions;
  Long_t i, j, k, size, sizeEtaBins, sizePhiBins;
  TBinMap::iterator itEtaBin;
  set< Double_t >::iterator itPhiBin;
  vector< Double_t > *phiBins;

  // read eta and phi bins
  param = GetParam("EtaPhiBins");
  size = param.GetSize();
  fBinMap.clear();
  fEtaBins.clear();
  fPhiBins.clear();
  for(i = 0; i < size/2; ++i)
  {
    paramEtaBins = param[i*2];
    sizeEtaBins = paramEtaBins.GetSize();
    paramPhiBins = param[i*2 + 1];
    sizePhiBins = paramPhiBins.GetSize();

    for(j = 0; j < sizeEtaBins; ++j)
    {
      for(k = 0; k < sizePhiBins; ++k)
      {
        fBinMap[paramEtaBins[j].GetDouble()].insert(paramPhiBins[k].GetDouble());
      }
    }
  }

  // for better performance we transform map of sets to parallel vectors:
  // vector< double > and vector< vector< double >* >
  for(itEtaBin = fBinMap.begin(); itEtaBin != fBinMap.end(); ++itEtaBin)
  {
    fEtaBins.push_back(itEtaBin->first);
    phiBins = new vector< double >(itEtaBin->second.size());
    fPhiBins.push_back(phiBins);
    phiBins->clear();
    for(itPhiBin = itEtaBin->second.begin(); itPhiBin != itEtaBin->second.end(); ++itPhiBin)
    {
      phiBins->push_back(*itPhiBin);
    }
  }

  // Eta x Phi smearing to be applied
  fEtaPhiRes = GetDouble("EtaPhiRes", 0.005);

  fTrackInputArray = ImportArray(GetString("TrackInputArray", "TrackMerger/tracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();

  fDenseChargedInputArray = ImportArray(GetString("DenseChargedInputArray", "DenseMergeTracks/tracks"));
  fItDenseChargedInputArray = fDenseChargedInputArray->MakeIterator();

  fTrackOutputArray = ExportArray(GetString("TrackOutputArray", "tracks"));
}

//------------------------------------------------------------------------------

void DenseTrackFilter::Finish()
{
  vector< vector< Double_t >* >::iterator itPhiBin;
  if(fItTrackInputArray) delete fItTrackInputArray;
  if(fItDenseChargedInputArray) delete fItDenseChargedInputArray;
  for(itPhiBin = fPhiBins.begin(); itPhiBin != fPhiBins.end(); ++itPhiBin)
  {
    delete *itPhiBin;
  }
}

//------------------------------------------------------------------------------

void DenseTrackFilter::Process()
{
  Candidate *candidate, *track, *bestTrack, *mother, *trackRef, *bestTrackRef;
  TLorentzVector position, momentum;
  Short_t etaBin, phiBin, flags;
  Int_t number, towerTrackHits;
  Long64_t towerHit, towerEtaPhi, hitEtaPhi;

  Double_t pt, ptmax, eta, phi;
  Bool_t matched;

  vector< Double_t >::iterator itEtaBin;
  vector< Double_t >::iterator itPhiBin;
  vector< Double_t > *phiBins;

  vector< Long64_t >::iterator itTowerHits;

  fTowerHits.clear();

  // loop over all tracks
  fItDenseChargedInputArray->Reset();
  number = -1;
  while((track = static_cast<Candidate*>(fItDenseChargedInputArray->Next())))
  {
    const TLorentzVector &trackPosition = track->Position;
    ++number;

    // find eta bin [1, fEtaBins.size - 1]
    itEtaBin = lower_bound(fEtaBins.begin(), fEtaBins.end(), trackPosition.Eta());
    if(itEtaBin == fEtaBins.begin() || itEtaBin == fEtaBins.end()) continue;
    etaBin = distance(fEtaBins.begin(), itEtaBin);

    // phi bins for given eta bin
    phiBins = fPhiBins[etaBin];

    // find phi bin [1, phiBins.size - 1]
    itPhiBin = lower_bound(phiBins->begin(), phiBins->end(), trackPosition.Phi());
    if(itPhiBin == phiBins->begin() || itPhiBin == phiBins->end()) continue;
    phiBin = distance(phiBins->begin(), itPhiBin);

    flags = 1;

    // make tower hit {16-bits for eta bin number, 16-bits for phi bin number, 8-bits for flags, 24-bits for track number}
    towerHit = (Long64_t(etaBin) << 48) | (Long64_t(phiBin) << 32) | (Long64_t(flags) << 24) | Long64_t(number);

    fTowerHits.push_back(towerHit);
  }

  // all hits are sorted first by eta bin number, then by phi bin number,
  // then by flags and then by particle or track number
  sort(fTowerHits.begin(), fTowerHits.end());

  // loop over all hits
  towerEtaPhi = 0;
  bestTrack = 0;
  fTower = 0;
  ptmax = 0.0;
  towerTrackHits = 0;

  for(itTowerHits = fTowerHits.begin(); itTowerHits != fTowerHits.end(); ++itTowerHits)
  {
    towerHit = (*itTowerHits);
    flags = (towerHit >> 24) & 0x00000000000000FFLL;
    number = (towerHit) & 0x0000000000FFFFFFLL;
    hitEtaPhi = towerHit >> 32;
    if(towerEtaPhi != hitEtaPhi)
    {
      // switch to next tower
      towerEtaPhi = hitEtaPhi;
      
      // saving track with highest pT that hit the tower
      if(towerTrackHits > 0)
      {
         bestTrackRef = static_cast<Candidate*>(bestTrack->GetCandidates()->At(0));
	 
         // find corresponding track in properly propagated tracks
         
	 
	 fItTrackInputArray->Reset();
	 matched = kFALSE;
	 int ntrack = 0;
	 while((track = static_cast<Candidate*>(fItTrackInputArray->Next())))
         {
           ntrack++;
           trackRef = static_cast<Candidate*>(track->GetCandidates()->At(0));
           if (trackRef->GetUniqueID() == bestTrackRef->GetUniqueID())
           {
             mother = track;
             candidate = static_cast<Candidate*>(track->Clone());
             pt = candidate->Momentum.Pt();
             eta = candidate->Momentum.Eta();
             phi = candidate->Momentum.Phi();
             eta = gRandom->Gaus(eta, fEtaPhiRes);
             phi = gRandom->Gaus(phi, fEtaPhiRes);
             candidate->Momentum.SetPtEtaPhiE(pt, eta, phi, pt*TMath::CosH(eta));
             candidate->AddCandidate(mother);
             fTrackOutputArray->Add(candidate);
             matched = kTRUE;
	   }
           if (matched) break;
	 }
         
         // find track 
      }

      ptmax = 0.0;
      towerTrackHits = 0;
      bestTrack = 0;
    }
    // check for track hits
    
    if(flags & 1)
    {
      ++towerTrackHits;
      track = static_cast<Candidate*>(fDenseChargedInputArray->At(number));
      momentum = track->Momentum;

      if (momentum.Pt() > ptmax) 
      {
        ptmax = momentum.Pt();
        bestTrack = track;
      }
      continue;
    }
  }
  
}
