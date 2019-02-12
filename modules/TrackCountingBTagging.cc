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

  fUse3D = GetBool("Use3D", false);

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

  Double_t jpx, jpy, jpz;
  Double_t dr, tpt;
  Double_t xd, yd, zd, d0, dd0, dz, ddz, sip;

  Int_t sign;

  Int_t count;

  // loop over all input jets
  fItJetInputArray->Reset();
  while((jet = static_cast<Candidate *>(fItJetInputArray->Next())))
  {
    const TLorentzVector &jetMomentum = jet->Momentum;
    jpx = jetMomentum.Px();
    jpy = jetMomentum.Py();
    jpz = jetMomentum.Pz();

    // loop over all input tracks
    fItTrackInputArray->Reset();
    count = 0;
    // stop once we have enough tracks
    while((track = static_cast<Candidate *>(fItTrackInputArray->Next())) and count < fNtracks)
    {
      const TLorentzVector &trkMomentum = track->Momentum;
      tpt = trkMomentum.Pt();
      if(tpt < fPtMin) continue;

      d0 = TMath::Abs(track->D0);
      if(d0 > fIPmax) continue;

      dr = jetMomentum.DeltaR(trkMomentum);
      if(dr > fDeltaR) continue;

      xd = track->Xd;
      yd = track->Yd;
      zd = track->Zd;
      dd0 = TMath::Abs(track->ErrorD0);
      dz = TMath::Abs(track->DZ);
      ddz = TMath::Abs(track->ErrorDZ);

      if(fUse3D)
      {
        sign = (jpx * xd + jpy * yd + jpz * zd > 0.0) ? 1 : -1;
        //add transverse and longitudinal significances in quadrature
        sip = sign * TMath::Sqrt(TMath::Power(d0 / dd0, 2) + TMath::Power(dz / ddz, 2));
      }
      else
      {
        sign = (jpx * xd + jpy * yd > 0.0) ? 1 : -1;
        sip = sign * d0 / TMath::Abs(dd0);
      }

      if(sip > fSigMin) count++;
    }

    // set BTag flag to true if count >= Ntracks
    jet->BTag |= (count >= fNtracks) << fBitNumber;
  }
}

//------------------------------------------------------------------------------
