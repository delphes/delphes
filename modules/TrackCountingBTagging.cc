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


/** \class TrackCountingBTagging
 *
 *  b-tagging algorithm based on counting tracks with large impact parameter
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "modules/TrackCountingBTagging.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

TrackCountingBTagging::TrackCountingBTagging() :
  fItTrackInputArray(0), fItJetInputArray(0)
{
}

//------------------------------------------------------------------------------

TrackCountingBTagging::~TrackCountingBTagging()
{
}

//------------------------------------------------------------------------------

void TrackCountingBTagging::Init()
{
  fBitNumber = GetInt("BitNumber", 0);

  fPtMin = GetDouble("TrackPtMin", 1.0);
  fDeltaR = GetDouble("DeltaR", 0.3);
  fIPmax = GetDouble("TrackIPMax", 2.0);

  fSigMin = GetDouble("SigMin", 6.5);
  fNtracks = GetInt("Ntracks", 3);

  // import input array(s)

  fTrackInputArray = ImportArray(GetString("TrackInputArray", "Calorimeter/eflowTracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();

  fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();
}

//------------------------------------------------------------------------------

void TrackCountingBTagging::Finish()
{
  if(fItTrackInputArray) delete fItTrackInputArray;
  if(fItJetInputArray) delete fItJetInputArray;
}

//------------------------------------------------------------------------------

void TrackCountingBTagging::Process()
{
  Candidate *jet, *track;

  Double_t jpx, jpy;
  Double_t dr, tpx, tpy, tpt;
  Double_t xd, yd, dxy, ddxy, ip, sip;

  Int_t sign;

  Int_t count;

  // loop over all input jets
  fItJetInputArray->Reset();
  while((jet = static_cast<Candidate*>(fItJetInputArray->Next())))
  {
    const TLorentzVector &jetMomentum = jet->Momentum;
    jpx = jetMomentum.Px();
    jpy = jetMomentum.Py();

    // loop over all input tracks
    fItTrackInputArray->Reset();
    count = 0;
    while((track = static_cast<Candidate*>(fItTrackInputArray->Next())))
    {
      const TLorentzVector &trkMomentum = track->Momentum;

      dr = jetMomentum.DeltaR(trkMomentum);

      tpt = trkMomentum.Pt();
      tpx = trkMomentum.Px();
      tpy = trkMomentum.Py();

      xd = track->Xd;
      yd = track->Yd;
      dxy = TMath::Abs(track->Dxy);
      ddxy = track->SDxy;

      if(tpt < fPtMin) continue;
      if(dr > fDeltaR) continue;
      if(dxy > fIPmax) continue;

      sign = (jpx*xd + jpy*yd > 0.0) ? 1 : -1;

      ip = sign*dxy;
      sip = ip / TMath::Abs(ddxy);

      if(sip > fSigMin) count++;
    }

    // set BTag flag to true if count >= Ntracks
    jet->BTag |= (count >= fNtracks) << fBitNumber;
  }
}

//------------------------------------------------------------------------------
